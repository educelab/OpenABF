/**
 * @example BenchmarkFlattening.cpp
 *
 * # Flattening benchmark
 *
 * Measures wall-clock runtime for three flattening configurations on one or
 * more mesh files and prints the results as a Markdown table:
 *
 * | Num. faces | ABF++ (s) | LSCM SparseLU (s) | LSCM LSCG (s) | HLSCM (s) |
 *
 * Usage:
 * @code
 *   openabf_example_benchmark mesh1.obj [mesh2.ply ...]
 * @endcode
 */
#include <chrono>
#include <filesystem>
#include <iomanip>
#include <iostream>

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
    if (argc < 2) {
        std::cerr << "Usage: " << fs::path(argv[0]).filename().string()
                  << " mesh1.(obj|ply) [mesh2 ...]\n";
        return EXIT_FAILURE;
    }

    // Solver type aliases
    using Mtx = Eigen::SparseMatrix<float>;
    using ABFMesh = OpenABF::detail::ABF::Mesh<float>;
    using ABF = OpenABF::ABFPlusPlus<float, ABFMesh>;
    using LU = Eigen::SparseLU<Mtx>;
    using LSCG = Eigen::LeastSquaresConjugateGradient<Mtx>;
    using LSCM_LU = OpenABF::AngleBasedLSCM<float, ABFMesh, LU>;
    using LSCM_LSCG = OpenABF::AngleBasedLSCM<float, ABFMesh, LSCG>;
    using HLSCM = OpenABF::HierarchicalLSCM<float, ABFMesh, LSCG>;

    // Table header
    std::cout << "| Mesh | Num. faces | ABF++ (s) | LSCM SparseLU (s) | LSCM LSCG (s) "
                 "| HLSCM (s) |\n";
    std::cout << "|------|-----------|-----------|-------------------|---------------|"
                 "----------|\n";

    for (int i = 1; i < argc; ++i) {
        const fs::path path = argv[i];

        // Load mesh once; clone for each solver run
        auto baseMesh = OpenABF::ReadMesh<ABFMesh>(path);
        const auto numFaces = baseMesh->num_faces();

        std::cerr << "Benchmarking " << path.filename().string() << " (" << numFaces
                  << " faces)...\n";

        // ABF++ — time the angle optimization only (shared across LSCM variants)
        double abfTime{0};
        {
            auto mesh = baseMesh->clone();
            abfTime = timeIt([&] {
                std::size_t iters{0};
                float grad{OpenABF::INF<float>};
                ABF::Compute(mesh, iters, grad);
            });
        }

        // Helper: run ABF++ then time a given LSCM variant
        auto runLSCM = [&](auto computeLSCM) -> double {
            auto mesh = baseMesh->clone();
            std::size_t iters{0};
            float grad{OpenABF::INF<float>};
            ABF::Compute(mesh, iters, grad);
            return timeIt([&] { computeLSCM(mesh); });
        };

        double luTime = runLSCM([](auto& m) { LSCM_LU::Compute(m); });
        double lscgTime = runLSCM([](auto& m) { LSCM_LSCG::Compute(m); });
        double hlscmTime = runLSCM([](auto& m) { HLSCM::Compute(m); });

        std::cout << std::fixed << std::setprecision(2);
        std::cout << "| " << path.filename().string() << " | " << numFaces << " | " << abfTime
                  << " | " << luTime << " | " << lscgTime << " | " << hlscmTime << " |\n";
        std::cout.flush();
    }

    return EXIT_SUCCESS;
}
