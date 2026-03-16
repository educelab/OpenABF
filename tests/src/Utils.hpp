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
/**
 * Construct a UV hemisphere (half-sphere open at the equator).
 *
 * Vertex 0 is the north pole. The last ring of vertices forms the equatorial
 * boundary. Produces rings * sectors + 1 vertices and
 * 2 * (rings - 1) * sectors + sectors faces.
 *
 * Example: ConstructHemisphere(12, 24) → 289 vertices, 552 faces.
 */
template <typename MeshType>
auto ConstructHemisphere(std::size_t rings, std::size_t sectors) -> typename MeshType::Pointer
{
    using T = float;
    auto mesh = MeshType::New();

    // North pole (vertex 0)
    mesh->insert_vertex(T(0), T(0), T(1));

    // Latitude rings from near-pole to equator
    for (std::size_t r = 1; r <= rings; ++r) {
        T phi = PI<T> / T(2) * T(r) / T(rings);  // 0 at pole, π/2 at equator
        T sinPhi = std::sin(phi);
        T cosPhi = std::cos(phi);
        for (std::size_t s = 0; s < sectors; ++s) {
            T theta = T(2) * PI<T> * T(s) / T(sectors);
            mesh->insert_vertex(sinPhi * std::cos(theta), sinPhi * std::sin(theta), cosPhi);
        }
    }

    // Faces: polar cap (triangles from pole to first ring)
    std::vector<std::vector<std::size_t>> allFaces;
    for (std::size_t s = 0; s < sectors; ++s) {
        std::size_t next = (s + 1) % sectors;
        allFaces.push_back({0, 1 + s, 1 + next});
    }

    // Faces: quad strips between consecutive rings
    for (std::size_t r = 0; r < rings - 1; ++r) {
        std::size_t base0 = 1 + r * sectors;
        std::size_t base1 = 1 + (r + 1) * sectors;
        for (std::size_t s = 0; s < sectors; ++s) {
            std::size_t next = (s + 1) % sectors;
            allFaces.push_back({base0 + s, base1 + s, base0 + next});
            allFaces.push_back({base0 + next, base1 + s, base1 + next});
        }
    }

    mesh->insert_faces(allFaces);
    return mesh;
}

/**
 * Construct a wavy surface: a grid with z-displacement
 *   z = 0.3 * sin(2π x / (cols-1)) * cos(2π y / (rows-1))
 *
 * Same connectivity as ConstructGrid but with spatially varying curvature.
 *
 * Example: ConstructWavySurface(20, 20) → 400 vertices, 722 faces.
 */
template <typename MeshType>
auto ConstructWavySurface(std::size_t rows, std::size_t cols) -> typename MeshType::Pointer
{
    using T = float;
    auto mesh = MeshType::New();

    for (std::size_t r = 0; r < rows; ++r) {
        for (std::size_t c = 0; c < cols; ++c) {
            T x = T(c);
            T y = T(r);
            T z = T(0.3) * std::sin(T(2) * PI<T> * x / T(cols - 1)) *
                  std::cos(T(2) * PI<T> * y / T(rows - 1));
            mesh->insert_vertex(x, y, z);
        }
    }

    std::vector<std::vector<std::size_t>> faces;
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

}  // namespace OpenABF::tests
