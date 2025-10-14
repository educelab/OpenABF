/**
 * @example CloneMesh.cpp
 *
 * # Mesh cloning demo
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
    const auto original = ABF::Mesh::New();
    original->insert_vertices(vertices);

    // Create the 3 pyramid triangles (the 4th is open)
    original->insert_faces({{1, 3, 0}, {3, 2, 0}, {3, 1, 2}});

    // Clone the original mesh
    auto flat = original->clone();

    // Run ABF++ on the clone
    std::size_t iters{0};
    float grad{OpenABF::INF<float>};
    ABF::Compute(flat, iters, grad);

    // Run LSCM
    LSCM::Compute(flat);

    // Write the meshes to verify the original is unaltered
    WriteMesh("openabf_example_clone_original.obj", original);
    WriteMesh("openabf_example_clone_flattened.obj", flat);
}