/**
 * @example BasicFlattening.cpp
 *
 * # Basic flattening demo
 *
 * Demonstrates flattening of a triangular pyramid mesh using ABF++ and LSCM.
 *
 * @note The `HalfEdgeMesh` class
 * [currently assumes](https://gitlab.com/educelab/OpenABF/-/issues/4) that the
 * surface has a boundary, is manifold, and that the winding order of all faces
 * is the same. Care should be taken that this assumption is not violated when
 * constructing your mesh.
 */
#include <iostream>

#include "OpenABF/OpenABF.hpp"

int main()
{
    // Define the flattening types
    using ABF = OpenABF::ABFPlusPlus<float>;
    using LSCM = OpenABF::AngleBasedLSCM<float, ABF::Mesh>;
    using Vec3f = OpenABF::Vec3f;

    // Pre-define vertices list
    const std::vector vertices = {
        Vec3f{0.f, 0.f, 0.f}, Vec3f{2.f, 0.f, 0.f},
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
        std::cout << v->idx << ": " << vertices[v->idx] << " -> " << v->pos
                  << std::endl;
    }

    // Write the flattened mesh
    WriteMesh("openabf_example_basic_flattening.obj", mesh);
}
