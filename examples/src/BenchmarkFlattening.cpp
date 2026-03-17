/**
 * @example BenchmarkFlattening.cpp
 *
 * # Flattening benchmark
 *
 * Measures wall-clock runtime for flattening configurations on one or more
 * mesh files and prints results as a Markdown table.  LSCM LSCG is timed at
 * 1 thread then at powers of 2 up to --threads N (default: hardware
 * concurrency), matching the format of volume-cartographer#123 with an added
 * HLSCM column.
 *
 * Usage:
 * @code
 *   openabf_example_benchmark [--threads N] mesh1.obj [mesh2.ply ...]
 * @endcode
 *
 * @note If Eigen was not compiled with OpenMP, all LSCG columns will report
 * the same time regardless of the requested thread count.
 */
#include <chrono>
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

template <class Fn>
auto timeIt(Fn&& fn) -> double
{
    auto t0 = Clock::now();
    fn();
    return Seconds(Clock::now() - t0).count();
}

auto main(const int argc, char* argv[]) -> int
{
    // Parse optional --threads N before mesh paths
    int meshArgStart = 1;
    int maxThreads = static_cast<int>(std::thread::hardware_concurrency());
    if (maxThreads < 1) {
        maxThreads = 1;
    }

    if (argc >= 3 && std::string(argv[1]) == "--threads") {
        maxThreads = std::stoi(argv[2]);
        meshArgStart = 3;
    }

    if (meshArgStart >= argc) {
        std::cerr << "Usage: " << fs::path(argv[0]).filename().string()
                  << " [--threads N] mesh1.(obj|ply) [mesh2 ...]\n";
        return EXIT_FAILURE;
    }

    // Build thread-count sequence: 1, 2, 4, 8, ... <= maxThreads
    std::vector<int> threadCounts;
    for (int t = 1; t <= maxThreads; t *= 2) {
        threadCounts.push_back(t);
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
    using ABFMesh = OpenABF::detail::ABF::Mesh<float>;
    using ABF = OpenABF::ABFPlusPlus<float, ABFMesh>;
    using LU = Eigen::SparseLU<Mtx>;
    using LSCG = Eigen::LeastSquaresConjugateGradient<Mtx>;
    using LSCM_LU = OpenABF::AngleBasedLSCM<float, ABFMesh, LU>;
    using LSCM_LSCG = OpenABF::AngleBasedLSCM<float, ABFMesh, LSCG>;
    using HLSCM = OpenABF::HierarchicalLSCM<float, ABFMesh, LSCG>;

    // Print table header
    std::cout << "| Mesh | Num. faces | ABF++ (s) | LSCM SparseLU (s)";
    for (int t : actualThreads) {
        std::cout << " | LSCM LSCG (" << t << "t) (s)";
    }
    std::cout << " | HLSCM (s) |\n";

    std::cout << "|------|-----------|-----------|------------------";
    for (int i = 0; i < static_cast<int>(actualThreads.size()); ++i) {
        std::cout << "-|------------------";
    }
    std::cout << "-|----------|\n";

    for (int i = meshArgStart; i < argc; ++i) {
        const fs::path path = argv[i];

        auto baseMesh = OpenABF::ReadMesh<ABFMesh>(path);
        const auto numFaces = baseMesh->num_faces();

        std::cerr << "Benchmarking " << path.filename().string() << " (" << numFaces
                  << " faces)...\n";

        // Helper: run ABF++ then time a given LSCM variant at a thread count
        auto runLSCM = [&](int threads, auto computeLSCM) -> double {
            auto mesh = baseMesh->clone();
            std::size_t iters{0};
            float grad{OpenABF::INF<float>};
            ABF::Compute(mesh, iters, grad);
            Eigen::setNbThreads(threads);
            double t = timeIt([&] { computeLSCM(mesh); });
            Eigen::setNbThreads(1);
            return t;
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

        double luTime = runLSCM(1, [](auto& m) { LSCM_LU::Compute(m); });

        std::vector<double> lscgTimes;
        for (int t : actualThreads) {
            lscgTimes.push_back(runLSCM(t, [](auto& m) { LSCM_LSCG::Compute(m); }));
        }

        double hlscmTime = runLSCM(1, [](auto& m) { HLSCM::Compute(m); });

        std::cout << std::fixed << std::setprecision(2);
        std::cout << "| " << path.filename().string() << " | " << numFaces << " | " << abfTime
                  << " | " << luTime;
        for (double lt : lscgTimes) {
            std::cout << " | " << lt;
        }
        std::cout << " | " << hlscmTime << " |\n";
        std::cout.flush();
    }

    return EXIT_SUCCESS;
}
