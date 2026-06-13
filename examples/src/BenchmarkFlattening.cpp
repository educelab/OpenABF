/**
 * @example BenchmarkFlattening.cpp
 *
 * # Flattening benchmark
 *
 * Measures wall-clock runtime for flattening configurations on one or more
 * meshes and prints results as a Markdown table.  LSCM CG and HLSCM variants
 * are timed at each specified thread count.
 *
 * Usage:
 * @code
 *   openabf_example_benchmark [OPTIONS] [mesh1.obj ...]
 * @endcode
 *
 * Options:
 *   --threads N [N ...]  Thread counts to benchmark (default: 1, 2, 4, …,
 *                        hardware concurrency). Accepts one or more values,
 *                        e.g. --threads 1 4 12.
 *   --output-dir DIR     Write one flattened OBJ per algorithm per mesh into
 *                        DIR.
 *   --builtin [MAX]      Benchmark the built-in wavy-surface sequence
 *                        (50k, 100k, 200k, 400k, 600k, 800k, 1M faces) up to
 *                        MAX faces (default: 1000000). Mesh files and
 *                        --builtin may be combined.
 *   --builtin-sphere [MAX]
 *                        Benchmark the built-in sphere-cap (hemisphere)
 *                        sequence using the same face-count progression as
 *                        --builtin. Sphere caps have constant positive
 *                        Gaussian curvature and exercise ABF+LSCM quality
 *                        in addition to performance.
 *
 * @note If Eigen was not compiled with OpenMP, all multi-thread columns will
 * report the same time as the 1-thread column.
 *
 * @note SparseLU can exhaust memory on large meshes. If a solver throws
 * std::bad_alloc or OpenABF::SolverException the cell is reported as "N/A"
 * and the run continues.
 */
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <thread>
#include <vector>

#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>

#include "OpenABF/OpenABF.hpp"

namespace fs = std::filesystem;
using Clock = std::chrono::steady_clock;
using Seconds = std::chrono::duration<double>;

using FloatT = double;
using ABFMesh = OpenABF::detail::ABF::Mesh<FloatT>;

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
    using T = FloatT;
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

/** Build a sphere cap with approximately targetFaces triangles.
 *
 * The cap spans colatitude [0, thetaMax] over the full longitude range;
 * thetaMax = π/2 yields a unit hemisphere. The cap has constant positive
 * Gaussian curvature, which forces ABF++ to take many iterations and
 * produces a harder LSCM system than the wavy surface.
 *
 * The pole is vertex 0; each successive ring has `sectors` vertices.
 * Produces sectors * (2*rings - 1) faces. Rings and sectors are chosen so
 * each triangle is roughly square: arc length per ring (thetaMax / rings)
 * ≈ arc length per sector at the boundary (2π * sin(thetaMax) / sectors).
 */
auto buildSphereCap(std::size_t targetFaces, FloatT thetaMax = OpenABF::PI<FloatT> / FloatT(2)) ->
    typename ABFMesh::Pointer
{
    using T = FloatT;
    // sectors / rings to keep triangles roughly square
    T aspect = T(2) * OpenABF::PI<T> * std::sin(thetaMax) / thetaMax;
    // sectors * (2*rings - 1) ≈ targetFaces, with sectors ≈ aspect * rings
    auto rings = static_cast<std::size_t>(
                     std::floor(std::sqrt(static_cast<double>(targetFaces) / (2.0 * aspect)))) +
                 1;
    auto sectors = static_cast<std::size_t>(std::floor(aspect * static_cast<double>(rings)));
    if (sectors < 3) {
        sectors = 3;
    }

    auto mesh = ABFMesh::New();
    mesh->insert_vertex(T(0), T(0), T(1));
    for (std::size_t r = 1; r <= rings; ++r) {
        T phi = thetaMax * T(r) / T(rings);
        T sinPhi = std::sin(phi);
        T cosPhi = std::cos(phi);
        for (std::size_t s = 0; s < sectors; ++s) {
            T theta = T(2) * OpenABF::PI<T> * T(s) / T(sectors);
            mesh->insert_vertex(sinPhi * std::cos(theta), sinPhi * std::sin(theta), cosPhi);
        }
    }
    std::vector<std::vector<std::size_t>> faces;
    faces.reserve(sectors + 2 * (rings - 1) * sectors);
    for (std::size_t s = 0; s < sectors; ++s) {
        std::size_t next = (s + 1) % sectors;
        faces.push_back({0, 1 + s, 1 + next});
    }
    for (std::size_t r = 0; r < rings - 1; ++r) {
        std::size_t base0 = 1 + r * sectors;
        std::size_t base1 = 1 + (r + 1) * sectors;
        for (std::size_t s = 0; s < sectors; ++s) {
            std::size_t next = (s + 1) % sectors;
            faces.push_back({base0 + s, base1 + s, base0 + next});
            faces.push_back({base0 + next, base1 + s, base1 + next});
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

/** Sentinel value indicating a solver failed (OOM or SolverException) */
static constexpr double kFailed = -1.0;

auto main(const int argc, char* argv[]) -> int
{
    fs::path outputDir;
    bool builtinEnabled = false;
    std::size_t builtinMax = 1'000'000;
    bool builtinSphereEnabled = false;
    std::size_t builtinSphereMax = 1'000'000;
    std::vector<fs::path> meshFiles;
    std::vector<int> threadCounts;

    for (int a = 1; a < argc; ++a) {
        std::string arg = argv[a];
        if (arg == "--threads") {
            // Consume one or more following numeric arguments
            while (a + 1 < argc) {
                std::string next = argv[a + 1];
                if (next.rfind("--", 0) == 0) {
                    break;
                }
                try {
                    threadCounts.push_back(std::stoi(next));
                    ++a;
                } catch (...) {
                    break;
                }
            }
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
        } else if (arg == "--builtin-sphere") {
            builtinSphereEnabled = true;
            if (a + 1 < argc) {
                try {
                    builtinSphereMax = std::stoull(argv[a + 1]);
                    ++a;
                } catch (...) {
                }
            }
        } else {
            meshFiles.emplace_back(argv[a]);
        }
    }

    if (!builtinEnabled && !builtinSphereEnabled && meshFiles.empty()) {
        std::cerr << "Usage: " << fs::path(argv[0]).filename().string()
                  << " [--threads N [N ...]] [--output-dir DIR]"
                     " [--builtin [MAX_FACES]] [--builtin-sphere [MAX_FACES]]"
                     " [mesh1.(obj|ply) ...]\n";
        return EXIT_FAILURE;
    }

    if (!outputDir.empty()) {
        fs::create_directories(outputDir);
    }

    // Build thread-count sequence: if explicit values were given use them;
    // otherwise default to 1, 2, 4, … <= hardware concurrency, plus the
    // hardware concurrency if not already a power of 2.
    if (threadCounts.empty()) {
        int maxThreads = static_cast<int>(std::thread::hardware_concurrency());
        if (maxThreads < 1) {
            maxThreads = 1;
        }
        for (int t = 1; t <= maxThreads; t *= 2) {
            threadCounts.push_back(t);
        }
        if (threadCounts.back() != maxThreads) {
            threadCounts.push_back(maxThreads);
        }
    }

    // Determine actual counts Eigen will use (clamped to 1 without OpenMP)
    std::vector<int> actualThreads;
    for (int t : threadCounts) {
        Eigen::setNbThreads(t);
        actualThreads.push_back(Eigen::nbThreads());
    }
    Eigen::setNbThreads(1);

    // Solver type aliases
    using Mtx = Eigen::SparseMatrix<FloatT>;
    using ABF = OpenABF::ABFPlusPlus<FloatT, ABFMesh>;
    using LU = Eigen::SparseLU<Mtx>;
    // Diagonal (Jacobi) preconditioner — Eigen's CG default.
    using CG_Diag = Eigen::ConjugateGradient<Mtx, Eigen::Lower | Eigen::Upper>;
    // IncompleteCholesky preconditioner — typically converges in fewer
    // iterations than Diagonal on the LSCM normal equations.
    using CG_IC = Eigen::ConjugateGradient<Mtx, Eigen::Lower | Eigen::Upper,
                                           Eigen::IncompleteCholesky<FloatT>>;
    using LSCM_LU = OpenABF::AngleBasedLSCM<FloatT, ABFMesh, LU>;
    using LSCM_CG_Diag = OpenABF::AngleBasedLSCM<FloatT, ABFMesh, CG_Diag>;
    using LSCM_CG_IC = OpenABF::AngleBasedLSCM<FloatT, ABFMesh, CG_IC>;
    using HLSCM_LSCG = OpenABF::HierarchicalLSCM<FloatT, ABFMesh>;
    using HLSCM_CG = OpenABF::HierarchicalLSCM<FloatT, ABFMesh, CG_Diag>;

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

    if (builtinSphereEnabled) {
        for (std::size_t n :
             {50'000UL, 100'000UL, 200'000UL, 400'000UL, 600'000UL, 800'000UL, 1'000'000UL}) {
            if (n > builtinSphereMax) {
                break;
            }
            auto mesh = buildSphereCap(n);
            auto label = "sphere~" + std::to_string(mesh->num_faces()) + "f";
            inputs.push_back({label, mesh});
        }
    }

    for (const auto& p : meshFiles) {
        inputs.push_back({p.filename().string(), OpenABF::ReadMesh<ABFMesh>(p)});
    }

    // Print table header
    std::cout << "| Mesh | Num. faces | ABF++ (s) | LSCM SparseLU (s)";
    for (int t : actualThreads) {
        std::cout << " | LSCM CG-Diag (" << t << "t) (s)";
    }
    for (int t : actualThreads) {
        std::cout << " | LSCM CG-IC (" << t << "t) (s)";
    }
    for (int t : actualThreads) {
        std::cout << " | HLSCM LSCG (" << t << "t) (s)";
    }
    for (int t : actualThreads) {
        std::cout << " | HLSCM CG (" << t << "t) (s)";
    }
    std::cout << " |\n";

    std::cout << "|------|-----------|-----------|------------------";
    for (std::size_t i = 0; i < 4 * actualThreads.size(); ++i) {
        std::cout << "-|------------------";
    }
    std::cout << "-|\n";

    // Helper to format a timed result (prints "N/A" for failures)
    auto fmtTime = [](double t) -> std::string {
        if (t < 0) {
            return "N/A";
        }
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(2) << t;
        return oss.str();
    };

    for (auto& input : inputs) {
        const auto& label = input.label;
        const auto& baseMesh = input.mesh;
        const auto numFaces = baseMesh->num_faces();
        std::cerr << "Benchmarking " << label << " (" << numFaces << " faces)...\n";

        // Helper: run ABF++ then time a given LSCM variant.
        // Returns {kFailed, nullptr} if the solver throws bad_alloc or
        // SolverException so the benchmark can print "N/A" and continue.
        auto runLSCM = [&](int threads,
                           auto computeLSCM) -> std::pair<double, typename ABFMesh::Pointer> {
            auto mesh = baseMesh->clone();
            std::size_t iters{0};
            FloatT grad{OpenABF::INF<FloatT>};
            ABF::Compute(mesh, iters, grad);
            Eigen::setNbThreads(threads);
            double t = kFailed;
            try {
                t = timeIt([&] { computeLSCM(mesh); });
            } catch (const std::bad_alloc&) {
                std::cerr << "  [OOM]\n";
            } catch (const OpenABF::SolverException& e) {
                std::cerr << "  [solver error: " << e.what() << "]\n";
            }
            Eigen::setNbThreads(1);
            if (t < 0) {
                return {kFailed, nullptr};
            }
            return {t, mesh};
        };

        // Warmup: one untimed run of each threaded solver to hot-load mesh
        // data into cache, preventing the 1-thread column from being penalized
        // by cold-start effects.
        runLSCM(1, [](auto& m) { LSCM_CG_Diag::Compute(m); });
        runLSCM(1, [](auto& m) { LSCM_CG_IC::Compute(m); });
        runLSCM(1, [](auto& m) { HLSCM_LSCG::Compute(m); });
        runLSCM(1, [](auto& m) { HLSCM_CG::Compute(m); });

        // ABF++ timing
        double abfTime{0};
        {
            auto mesh = baseMesh->clone();
            abfTime = timeIt([&] {
                std::size_t iters{0};
                FloatT grad{OpenABF::INF<FloatT>};
                ABF::Compute(mesh, iters, grad);
            });
        }

        auto [luTime, luMesh] = runLSCM(1, [](auto& m) { LSCM_LU::Compute(m); });

        std::vector<double> cgDiagTimes;
        typename ABFMesh::Pointer cgDiagMesh;
        for (int t : actualThreads) {
            auto [time, mesh] = runLSCM(t, [](auto& m) { LSCM_CG_Diag::Compute(m); });
            cgDiagTimes.push_back(time);
            if (t == actualThreads.front()) {
                cgDiagMesh = mesh;
            }
        }

        std::vector<double> cgIcTimes;
        typename ABFMesh::Pointer cgIcMesh;
        for (int t : actualThreads) {
            auto [time, mesh] = runLSCM(t, [](auto& m) { LSCM_CG_IC::Compute(m); });
            cgIcTimes.push_back(time);
            if (t == actualThreads.front()) {
                cgIcMesh = mesh;
            }
        }

        std::vector<double> hlscmLscgTimes;
        typename ABFMesh::Pointer hlscmLscgMesh;
        for (int t : actualThreads) {
            auto [time, mesh] = runLSCM(t, [](auto& m) { HLSCM_LSCG::Compute(m); });
            hlscmLscgTimes.push_back(time);
            if (t == actualThreads.front()) {
                hlscmLscgMesh = mesh;
            }
        }

        std::vector<double> hlscmCgTimes;
        typename ABFMesh::Pointer hlscmCgMesh;
        for (int t : actualThreads) {
            auto [time, mesh] = runLSCM(t, [](auto& m) { HLSCM_CG::Compute(m); });
            hlscmCgTimes.push_back(time);
            if (t == actualThreads.front()) {
                hlscmCgMesh = mesh;
            }
        }

        // Optionally write flattened meshes
        if (!outputDir.empty()) {
            // Sanitize label for use as filename stem
            auto stem = label;
            for (auto& ch : stem) {
                if (ch == '~' || ch == ' ') {
                    ch = '_';
                }
            }
            OpenABF::WriteMesh(outputDir / (stem + "_3d.obj"), baseMesh);
            if (luMesh) {
                OpenABF::WriteMesh(outputDir / (stem + "_lscm_lu.obj"), luMesh);
            }
            if (cgDiagMesh) {
                OpenABF::WriteMesh(outputDir / (stem + "_lscm_cg_diag.obj"), cgDiagMesh);
            }
            if (cgIcMesh) {
                OpenABF::WriteMesh(outputDir / (stem + "_lscm_cg_ic.obj"), cgIcMesh);
            }
            if (hlscmLscgMesh) {
                OpenABF::WriteMesh(outputDir / (stem + "_hlscm_lscg.obj"), hlscmLscgMesh);
            }
            if (hlscmCgMesh) {
                OpenABF::WriteMesh(outputDir / (stem + "_hlscm_cg.obj"), hlscmCgMesh);
            }
        }

        std::cout << "| " << label << " | " << numFaces << " | " << std::fixed
                  << std::setprecision(2) << abfTime << " | " << fmtTime(luTime);
        for (double ct : cgDiagTimes) {
            std::cout << " | " << fmtTime(ct);
        }
        for (double ct : cgIcTimes) {
            std::cout << " | " << fmtTime(ct);
        }
        for (double ht : hlscmLscgTimes) {
            std::cout << " | " << fmtTime(ht);
        }
        for (double ht : hlscmCgTimes) {
            std::cout << " | " << fmtTime(ht);
        }
        std::cout << " |\n";
        std::cout.flush();
    }

    return EXIT_SUCCESS;
}
