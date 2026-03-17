/**
 * @example BenchmarkFlattening.cpp
 *
 * # Flattening benchmark
 *
 * Measures wall-clock runtime for flattening configurations on one or more
 * meshes and prints results as a Markdown table.  LSCM CG and HLSCM CG are
 * timed at 1 thread then at powers of 2 up to --threads N.
 *
 * Usage:
 * @code
 *   openabf_example_benchmark [OPTIONS] [mesh1.obj ...]
 * @endcode
 *
 * Options:
 *   --threads N        Max thread count for LSCG/HLSCM columns (default:
 *                      hardware concurrency). Columns run at 1, 2, 4, … N.
 *   --output-dir DIR   Write one flattened OBJ per algorithm per mesh into DIR.
 *   --builtin [MAX]    Benchmark the built-in wavy-surface sequence
 *                      (50k, 100k, 200k, 400k, 600k, 800k, 1M faces) up to MAX
 *                      faces (default: 1000000). Mesh files and --builtin may
 *                      be combined.
 *
 * @note If Eigen was not compiled with OpenMP, all multi-thread columns will
 * report the same time as the 1-thread column.
 */
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <thread>
#include <vector>

#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>

#include "OpenABF/OpenABF.hpp"

namespace fs = std::filesystem;
using Clock = std::chrono::steady_clock;
using Seconds = std::chrono::duration<double>;

using ABFMesh = OpenABF::detail::ABF::Mesh<float>;

template <class Fn>
auto timeIt(Fn&& fn) -> double
{
    auto t0 = Clock::now();
    fn();
    return Seconds(Clock::now() - t0).count();
}

/** Build a wavy surface with approximately targetFaces triangles.
 *
 * Computes a square-ish vertex grid with rows = cols = floor(sqrt(N/2)) + 1,
 * giving 2*(rows-1)*(cols-1) actual triangles.
 */
auto buildWavySurface(std::size_t targetFaces) -> typename ABFMesh::Pointer
{
    using T = float;
    auto n = static_cast<std::size_t>(std::floor(std::sqrt(targetFaces / 2.0))) + 1;
    std::size_t rows = n, cols = n;

    auto mesh = ABFMesh::New();
    for (std::size_t r = 0; r < rows; ++r) {
        for (std::size_t c = 0; c < cols; ++c) {
            T x = T(c);
            T y = T(r);
            T z = T(0.3) * std::sin(T(2) * OpenABF::PI<T> * x / T(cols - 1)) *
                  std::cos(T(2) * OpenABF::PI<T> * y / T(rows - 1));
            mesh->insert_vertex(x, y, z);
        }
    }
    std::vector<std::vector<std::size_t>> faces;
    faces.reserve(2 * (rows - 1) * (cols - 1));
    for (std::size_t r = 0; r < rows - 1; ++r) {
        for (std::size_t c = 0; c < cols - 1; ++c) {
            auto v0 = r * cols + c;
            auto v1 = r * cols + c + 1;
            auto v2 = (r + 1) * cols + c;
            auto v3 = (r + 1) * cols + c + 1;
            faces.push_back({v0, v2, v1});
            faces.push_back({v1, v2, v3});
        }
    }
    mesh->insert_faces(faces);
    return mesh;
}

/** Benchmark inputs: either a file path or a synthetic face count */
struct BenchInput {
    std::string label;
    typename ABFMesh::Pointer mesh;
};

auto main(const int argc, char* argv[]) -> int
{
    int maxThreads = static_cast<int>(std::thread::hardware_concurrency());
    if (maxThreads < 1) {
        maxThreads = 1;
    }
    fs::path outputDir;
    bool builtinEnabled = false;
    std::size_t builtinMax = 1'000'000;
    std::vector<fs::path> meshFiles;

    for (int a = 1; a < argc; ++a) {
        std::string arg = argv[a];
        if (arg == "--threads" && a + 1 < argc) {
            maxThreads = std::stoi(argv[++a]);
        } else if (arg == "--output-dir" && a + 1 < argc) {
            outputDir = argv[++a];
        } else if (arg == "--builtin") {
            builtinEnabled = true;
            // Optional next arg: max face count (if it parses as a number)
            if (a + 1 < argc) {
                try {
                    builtinMax = std::stoull(argv[a + 1]);
                    ++a;
                } catch (...) {
                }
            }
        } else {
            meshFiles.emplace_back(argv[a]);
        }
    }

    if (!builtinEnabled && meshFiles.empty()) {
        std::cerr << "Usage: " << fs::path(argv[0]).filename().string()
                  << " [--threads N] [--output-dir DIR] [--builtin [MAX_FACES]]"
                     " [mesh1.(obj|ply) ...]\n";
        return EXIT_FAILURE;
    }

    if (!outputDir.empty()) {
        fs::create_directories(outputDir);
    }

    // Build thread-count sequence: 1, 2, 4, … <= maxThreads, plus maxThreads
    // itself if not already a power of 2.
    std::vector<int> threadCounts;
    for (int t = 1; t <= maxThreads; t *= 2) {
        threadCounts.push_back(t);
    }
    if (threadCounts.back() != maxThreads) {
        threadCounts.push_back(maxThreads);
    }

    // Determine actual counts Eigen will use (clamped to 1 without OpenMP)
    std::vector<int> actualThreads;
    for (int t : threadCounts) {
        Eigen::setNbThreads(t);
        actualThreads.push_back(Eigen::nbThreads());
    }
    Eigen::setNbThreads(1);

    // Solver type aliases
    using Mtx = Eigen::SparseMatrix<float>;
    using ABF = OpenABF::ABFPlusPlus<float, ABFMesh>;
    using LU = Eigen::SparseLU<Mtx>;
    using CG = Eigen::ConjugateGradient<Mtx>;
    using LSCM_LU = OpenABF::AngleBasedLSCM<float, ABFMesh, LU>;
    using LSCM_CG = OpenABF::AngleBasedLSCM<float, ABFMesh, CG>;
    using HLSCM = OpenABF::HierarchicalLSCM<float, ABFMesh, CG>;

    // Assemble benchmark inputs
    std::vector<BenchInput> inputs;

    if (builtinEnabled) {
        // Standard sequence matching volume-cartographer#123
        for (std::size_t n :
             {50'000UL, 100'000UL, 200'000UL, 400'000UL, 600'000UL, 800'000UL, 1'000'000UL}) {
            if (n > builtinMax) {
                break;
            }
            auto mesh = buildWavySurface(n);
            auto label = "wavy~" + std::to_string(mesh->num_faces()) + "f";
            inputs.push_back({label, mesh});
        }
    }

    for (const auto& p : meshFiles) {
        inputs.push_back({p.filename().string(), OpenABF::ReadMesh<ABFMesh>(p)});
    }

    // Print table header
    std::cout << "| Mesh | Num. faces | ABF++ (s) | LSCM SparseLU (s)";
    for (int t : actualThreads) {
        std::cout << " | LSCM CG (" << t << "t) (s)";
    }
    for (int t : actualThreads) {
        std::cout << " | HLSCM CG (" << t << "t) (s)";
    }
    std::cout << " |\n";

    std::cout << "|------|-----------|-----------|------------------";
    for (std::size_t i = 0; i < 2 * actualThreads.size(); ++i) {
        std::cout << "-|------------------";
    }
    std::cout << "-|\n";

    for (auto& input : inputs) {
        const auto& label = input.label;
        const auto& baseMesh = input.mesh;
        const auto numFaces = baseMesh->num_faces();
        std::cerr << "Benchmarking " << label << " (" << numFaces << " faces)...\n";

        // Helper: run ABF++ then time a given LSCM variant; return {time, mesh}
        auto runLSCM = [&](int threads,
                           auto computeLSCM) -> std::pair<double, typename ABFMesh::Pointer> {
            auto mesh = baseMesh->clone();
            std::size_t iters{0};
            float grad{OpenABF::INF<float>};
            ABF::Compute(mesh, iters, grad);
            Eigen::setNbThreads(threads);
            double t = timeIt([&] { computeLSCM(mesh); });
            Eigen::setNbThreads(1);
            return {t, mesh};
        };

        // ABF++ timing
        double abfTime{0};
        {
            auto mesh = baseMesh->clone();
            abfTime = timeIt([&] {
                std::size_t iters{0};
                float grad{OpenABF::INF<float>};
                ABF::Compute(mesh, iters, grad);
            });
        }

        auto [luTime, luMesh] = runLSCM(1, [](auto& m) { LSCM_LU::Compute(m); });

        std::vector<double> cgTimes;
        typename ABFMesh::Pointer cgMesh;
        for (int t : actualThreads) {
            auto [time, mesh] = runLSCM(t, [](auto& m) { LSCM_CG::Compute(m); });
            cgTimes.push_back(time);
            if (t == 1) {
                cgMesh = mesh;
            }
        }

        std::vector<double> hlscmTimes;
        typename ABFMesh::Pointer hlscmMesh;
        for (int t : actualThreads) {
            auto [time, mesh] = runLSCM(t, [](auto& m) { HLSCM::Compute(m); });
            hlscmTimes.push_back(time);
            if (t == 1) {
                hlscmMesh = mesh;
            }
        }

        // Optionally write flattened meshes (1-thread; threading doesn't affect output)
        if (!outputDir.empty()) {
            // Sanitize label for use as filename stem
            auto stem = label;
            for (auto& ch : stem) {
                if (ch == '~' || ch == ' ') {
                    ch = '_';
                }
            }
            OpenABF::WriteMesh(outputDir / (stem + "_lscm_lu.obj"), luMesh);
            OpenABF::WriteMesh(outputDir / (stem + "_lscm_cg.obj"), cgMesh);
            OpenABF::WriteMesh(outputDir / (stem + "_hlscm_cg.obj"), hlscmMesh);
        }

        std::cout << std::fixed << std::setprecision(2);
        std::cout << "| " << label << " | " << numFaces << " | " << abfTime << " | " << luTime;
        for (double ct : cgTimes) {
            std::cout << " | " << ct;
        }
        for (double ht : hlscmTimes) {
            std::cout << " | " << ht;
        }
        std::cout << " |\n";
        std::cout.flush();
    }

    return EXIT_SUCCESS;
}
