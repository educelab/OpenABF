#pragma once

#include "OpenABF/OpenABF.hpp"

namespace OpenABF::tests
{
/** Construct a triangular pyramid with an open bottom */
template <typename MeshType>
auto ConstructPyramid() -> typename MeshType::Pointer
{
    auto mesh = MeshType::New();
    mesh->insert_vertex(0, 0, 0);
    mesh->insert_vertex(2, 0, 0);
    mesh->insert_vertex(1, std::sqrt(3), 0);
    mesh->insert_vertex(1, std::sqrt(3) / 3, std::sqrt(6) * 2 / 3);

    mesh->insert_faces({{1, 3, 0}, {3, 2, 0}, {3, 1, 2}});
    return mesh;
}

/**
 * Construct a flat rectangular grid mesh, triangulated by splitting each quad
 * along its diagonal.  A (rows x cols) vertex grid produces (rows-1)*(cols-1)*2
 * triangles.  Interior vertex count = (rows-2)*(cols-2).
 *
 * Example: ConstructGrid<MeshType>(4, 4) → 16 vertices, 18 faces, 4 interior vertices.
 */
template <typename MeshType>
auto ConstructGrid(std::size_t rows, std::size_t cols) -> typename MeshType::Pointer
{
    auto mesh = MeshType::New();
    for (std::size_t r = 0; r < rows; ++r) {
        for (std::size_t c = 0; c < cols; ++c) {
            mesh->insert_vertex(static_cast<float>(c), static_cast<float>(r), 0.f);
        }
    }
    for (std::size_t r = 0; r < rows - 1; ++r) {
        for (std::size_t c = 0; c < cols - 1; ++c) {
            auto v0 = r * cols + c;
            auto v1 = r * cols + c + 1;
            auto v2 = (r + 1) * cols + c;
            auto v3 = (r + 1) * cols + c + 1;
            mesh->insert_faces({{v0, v2, v1}, {v1, v2, v3}});
        }
    }
    return mesh;
}
}  // namespace OpenABF::tests
