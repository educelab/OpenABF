#include <gtest/gtest.h>

#include "OpenABF/OpenABF.hpp"

using namespace OpenABF;

template <typename T>
auto ConstructPyramid() -> typename HalfEdgeMesh<T>::Pointer
{
    auto mesh = HalfEdgeMesh<T>::New();
    mesh->insert_vertex(0, 0, 0);
    mesh->insert_vertex(2, 0, 0);
    mesh->insert_vertex(1, std::sqrt(3), 0);
    mesh->insert_vertex(1, std::sqrt(3) / 3, std::sqrt(6) * 2 / 3);

    mesh->insert_face(0, 3, 1);
    mesh->insert_face(0, 2, 3);
    mesh->insert_face(2, 1, 3);
    return mesh;
}

TEST(HalfEdgeMesh, BuildMesh) { EXPECT_NO_THROW(ConstructPyramid<float>()); }

TEST(HalfEdgeMesh, FaceVertexOutOfBounds)
{
    auto mesh = HalfEdgeMesh<float>::New();
    EXPECT_THROW(mesh->insert_face(0, 1, 2), std::out_of_range);
}

TEST(HalfEdgeMesh, InsertNonManifoldEdge)
{
    auto mesh = ConstructPyramid<float>();
    EXPECT_NO_THROW(
        mesh->insert_vertex(1, std::sqrt(3) / 3, -std::sqrt(6) * 2 / 3));
    EXPECT_THROW(mesh->insert_face(0, 3, 4), MeshException);
}

TEST(HalfEdgeMesh, HasBoundary)
{
    // Construct open pyramid
    auto mesh = ConstructPyramid<float>();
    EXPECT_TRUE(HasBoundary(mesh));

    // Close pyramid
    mesh->insert_face(2, 0, 1);
    EXPECT_FALSE(HasBoundary(mesh));
}

TEST(HalfEdgeMesh, IsManifold)
{
    // Construct open pyramid
    auto mesh = ConstructPyramid<float>();
    EXPECT_TRUE(IsManifold(mesh));

    // Close pyramid
    mesh->insert_face(2, 0, 1);
    EXPECT_TRUE(IsManifold(mesh));

    // Complicated geometry tests
    // Simple triangle with hole
    // Not manifold because hole touches outer border
    mesh = HalfEdgeMesh<float>::New();
    mesh->insert_vertex(0, 0, 0);
    mesh->insert_vertex(1, 0, 0);
    mesh->insert_vertex(2, 0, 0);
    mesh->insert_vertex(0.5F, std::sqrt(3) / 2, 0);
    mesh->insert_vertex(1.5F, std::sqrt(3) / 2, 0);
    mesh->insert_vertex(1, std::sqrt(3), 0);
    mesh->insert_face(0, 3, 1);
    mesh->insert_face(1, 4, 2);
    mesh->insert_face(3, 5, 4);
    EXPECT_FALSE(IsManifold(mesh));

    // Add triangles to make hole not border-adjacent, is manifold
    mesh->insert_vertex(-0.5, 1, 0);
    mesh->insert_vertex(0.5, -std::sqrt(3), 0);
    mesh->insert_vertex(2.5, 1, 0);
    mesh->insert_face(0, 6, 3);
    mesh->insert_face(3, 6, 5);
    mesh->insert_face(4, 5, 7);
    mesh->insert_face(4, 7, 2);
    mesh->insert_face(1, 2, 8);
    mesh->insert_face(0, 1, 8);
    EXPECT_TRUE(IsManifold(mesh));

    // Add some non-manifold geometry that can't be detected by edge pairing
    mesh->insert_vertex(1, std::sqrt(3) / 3, std::sqrt(6) * 2 / 3);
    EXPECT_NO_THROW(mesh->insert_face(1, 3, 9));
    EXPECT_NO_THROW(mesh->insert_face(6, 9, 5));
    EXPECT_FALSE(IsManifold(mesh));
}