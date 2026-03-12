/**
 * @example AlternativeSolvers.cpp
 *
 * # Alternative solvers demo
 *
 * Demonstrates flattening using a conjugate gradient solver instead of
 * SparseLU.
 */
#include <iostream>

#include <Eigen/IterativeLinearSolvers>

#include "OpenABF/OpenABF.hpp"

auto main() -> int
{
    // Define the matrix, mesh, and vec classes
    using Mtx = Eigen::SparseMatrix<float>;
    using Mesh = OpenABF::detail::ABF::Mesh<float>;
    using Vec3f = OpenABF::Vec3f;
    // Define the solvers and methods
    using ABFSolver = Eigen::ConjugateGradient<Mtx, Eigen::Lower | Eigen::Upper>;
    using ABF = OpenABF::ABFPlusPlus<float, Mesh, ABFSolver>;
    using LSCMSolver = Eigen::LeastSquaresConjugateGradient<Mtx>;
    using LSCM = OpenABF::AngleBasedLSCM<float, ABF::Mesh, LSCMSolver>;

    // Pre-define vertices list
    const std::vector vertices = {Vec3f{0.f, 0.f, 0.f}, Vec3f{2.f, 0.f, 0.f},
                                  Vec3f{1.f, std::sqrt(3.f), 0.f},
                                  Vec3f{1.f, std::sqrt(3.f) / 3.f, std::sqrt(6.f) * 2.f / 3.f}};

    // Create the 4 pyramid vertices
    auto mesh = ABF::Mesh::New();
    mesh->insert_vertices(vertices);

    // Create the 3 pyramid triangles (the 4th is open)
    mesh->insert_faces({{1, 3, 0}, {3, 2, 0}, {3, 1, 2}});

    // Run ABF++
    std::size_t iters{0};
    float grad{OpenABF::INF<float>};
    ABF::Compute(mesh, iters, grad);
    std::cout << "ABF++ Final gradient: " << grad << std::endl;
    std::cout << "ABF++ Iterations: " << iters << std::endl;

    // Run LSCM
    LSCM::Compute(mesh);

    // Print flattened positions
    for (const auto& v : mesh->vertices()) {
        std::cout << v->idx << ": " << vertices[v->idx] << " -> " << v->pos << std::endl;
    }

    // Write the flattened mesh
    WriteMesh("openabf_example_alternative_solvers.obj", mesh);
}
