/**
 * @example SplitPath.cpp
 *
 * # Path splitting demo
 *
 * Demonstrates flattening of a closed, triangular pyramid by introducing a
 * border/seam along two edges.
 *
 * @see OpenABF::HalfEdgeMesh::split_path
 */
#include <iomanip>
#include <iostream>

#include "OpenABF/OpenABF.hpp"

int main()
{
    // Define the flattening types
    using ABF = OpenABF::ABFPlusPlus<float>;
    using LSCM = OpenABF::AngleBasedLSCM<float, ABF::Mesh>;

    // Create the 4 pyramid vertices
    auto mesh = ABF::Mesh::New();
    mesh->insert_vertex(0, 0, 0);
    mesh->insert_vertex(2, 0, 0);
    mesh->insert_vertex(1, std::sqrt(3), 0);
    mesh->insert_vertex(1, std::sqrt(3) / 3, std::sqrt(6) * 2 / 3);

    // Create the 4 pyramid triangles
    mesh->insert_faces({{1, 3, 0}, {3, 2, 0}, {3, 1, 2}, {1, 0, 2}});

    // Check manifoldness and boundary condition
    std::cout << std::boolalpha;
    std::cout << "Is manifold: " << IsManifold(mesh) << std::endl;
    std::cout << "Has boundary: " << HasBoundary(mesh) << std::endl;
    std::cout << "Edges: " << std::endl;
    for (const auto& e : mesh->edges()) {
        std::cout << "  " << e->vertex->idx << " -> " << e->pair->vertex->idx
                  << std::endl;
    }

    // Insert a border seam between vertices 0 -> 1 -> 2
    std::cout << "Splitting path 0 -> 1 -> 2" << std::endl;
    mesh->split_path({0, 1, 2});

    // Reprint manifoldness and boundary condition
    std::cout << "Is manifold: " << IsManifold(mesh) << std::endl;
    std::cout << "Has boundary: " << HasBoundary(mesh) << std::endl;
    std::cout << "Edges: " << std::endl;
    for (const auto& e : mesh->edges()) {
        std::cout << "  " << e->vertex->idx << " -> " << e->pair->vertex->idx
                  << std::endl;
    }

    // Run ABF++
    std::size_t iters{0};
    float grad{OpenABF::INF<float>};
    ABF::Compute(mesh, iters, grad);
    std::cout << "ABF++ Final gradient: " << grad << std::endl;
    std::cout << "ABF++ Iterations: " << iters << std::endl;

    // Run LSCM
    LSCM::Compute(mesh);

    // Write the flattened mesh
    WriteMesh("openabf_example_split_path.obj", mesh);
}
