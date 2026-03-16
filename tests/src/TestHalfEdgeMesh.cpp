#include <gtest/gtest.h>

#include "OpenABF/OpenABF.hpp"
#include "Utils.hpp"

using namespace OpenABF;
using namespace OpenABF::tests;

using MeshType = HalfEdgeMesh<float>;

TEST(HalfEdgeMesh, BuildMesh) { EXPECT_NO_THROW(ConstructPyramid<MeshType>()); }

TEST(HalfEdgeMesh, FaceVertexOutOfBounds)
{
    const auto mesh = MeshType::New();
    EXPECT_THROW(mesh->insert_face(0, 1, 2), std::out_of_range);
}

TEST(HalfEdgeMesh, IterateBoundary)
{
    // Single face mesh
    const auto mesh = MeshType::New();
    mesh->insert_vertex(0, 0, 0);
    mesh->insert_vertex(1, 0, 0);
    mesh->insert_vertex(0, 1, 0);
    mesh->insert_vertex(1, 1, 0);
    mesh->insert_faces({{0, 2, 1}});

    std::vector<std::size_t> indices;
    std::vector<std::size_t> expected{5, 3, 1};
    auto boundary = mesh->boundaries()[0];
    for (const auto& e : boundary) {
        indices.push_back(e->idx);
    }
    EXPECT_EQ(indices, expected);

    // Double face mesh
    mesh->insert_faces({{1, 2, 3}});
    indices.clear();
    expected = {5, 9, 7, 1};
    boundary = mesh->boundaries()[0];
    for (const auto& e : boundary) {
        indices.push_back(e->idx);
    }
    EXPECT_EQ(indices, expected);
}

TEST(HalfEdgeMesh, ConnectedComponents)
{
    // Mesh with a single CC
    const auto mesh = MeshType::New();
    mesh->insert_vertex(0, 0, 0);
    mesh->insert_vertex(1, 0, 0);
    mesh->insert_vertex(0, 1, 0);
    mesh->insert_faces({{0, 2, 1}});
    auto ccs = mesh->connected_components();
    EXPECT_EQ(ccs.size(), 1);
    EXPECT_EQ(ccs.size(), mesh->num_connected_components());

    // Add a second CC
    mesh->insert_vertex(1, 1, 1);
    mesh->insert_vertex(0, 1, 1);
    mesh->insert_vertex(1, 0, 1);
    mesh->insert_faces({{3, 5, 4}});
    ccs = mesh->connected_components();
    EXPECT_EQ(ccs.size(), 2);
    EXPECT_EQ(ccs.size(), mesh->num_connected_components());
}

TEST(HalfEdgeMesh, CheckWindingOrder)
{
    // Build basic mesh
    const auto mesh = MeshType::New();
    mesh->insert_vertex(0, 0, 0);
    mesh->insert_vertex(1, 0, 0);
    mesh->insert_vertex(0, 1, 0);
    mesh->insert_vertex(1, 1, 0);
    mesh->insert_faces({{0, 2, 1}});

    // Add adjacent mesh with wrong winding order
    std::size_t fid{0};
    EXPECT_NO_THROW(fid = mesh->insert_faces({{2, 1, 3}})[0]);
    std::vector<std::size_t> vids;
    std::vector<std::size_t> expected{1, 2, 3};
    for (const auto& e : *mesh->faces()[fid]) {
        vids.emplace_back(e->vertex->idx);
    }
    EXPECT_EQ(vids, expected);

    // Add face with bad winding order that only shares vertex
    mesh->insert_vertex(1, 0, 1);
    mesh->insert_vertex(1, 1, 1);
    EXPECT_THROW(mesh->insert_faces({{1, 5, 4}}), MeshException);

    // Add face between the two winding orders (currently not handled)
    EXPECT_THROW(mesh->insert_faces({{1, 3, 4}}), MeshException);
}

TEST(HalfEdgeMesh, InsertNonManifoldEdge)
{
    // Build a manifold mesh
    const auto mesh = MeshType::New();
    mesh->insert_vertex(0, 0, 0);
    mesh->insert_vertex(1, 0, 0);
    mesh->insert_vertex(0, 1, 0);
    mesh->insert_vertex(1, 1, 0);
    mesh->insert_vertex(1, 1, 1);

    mesh->insert_faces({{0, 2, 1}, {1, 2, 3}});

    // Try to add a non-manifold edge
    EXPECT_THROW(mesh->insert_faces({{1, 2, 4}}), MeshException);
}

TEST(HalfEdgeMesh, HasBoundary)
{
    // Construct open pyramid
    const auto mesh = ConstructPyramid<MeshType>();
    EXPECT_TRUE(HasBoundary(mesh));

    // Close pyramid
    mesh->insert_face(2, 0, 1);
    EXPECT_FALSE(HasBoundary(mesh));
}

TEST(HalfEdgeMesh, IsManifold)
{
    // Construct open pyramid
    auto mesh = ConstructPyramid<MeshType>();
    EXPECT_TRUE(IsManifold(mesh));

    // Close pyramid
    mesh->insert_faces({{2, 0, 1}});
    EXPECT_TRUE(IsManifold(mesh));

    // NOTE: The following manifold problems are caught when adding faces
    // to the mesh using add_faces. Leaving these tests here to test whether
    // we can detect some non-manifold constructions that are generally avoided
    // (e.g. creating a hole which touches the outer boundary).

    // Complicated geometry tests
    // Simple triangle with hole
    // Not manifold because hole touches outer border
    //      /\
    //     /  \
    //    /____\
    //   /\    /\
    //  /  \  /  \
    // /____\/____\

    mesh = MeshType::New();
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

TEST(HalfEdgeMesh, SimpleSplit)
{
    // Build a square mesh
    const auto mesh = MeshType::New();
    mesh->insert_vertex(0, 0, 0);
    mesh->insert_vertex(1, 0, 0);
    mesh->insert_vertex(0, 1, 0);
    mesh->insert_vertex(1, 1, 0);
    EXPECT_EQ(mesh->num_connected_components(), 0);

    const auto fids = mesh->insert_faces({{0, 2, 1}, {1, 2, 3}});
    EXPECT_EQ(mesh->num_connected_components(), 1);

    // Split along the diagonal
    mesh->split_path({1, 2});
    EXPECT_EQ(mesh->num_connected_components(), 2);
}

TEST(HalfEdgeMesh, TwoStageSplit)
{
    // Build a mesh
    const auto mesh = MeshType::New();
    mesh->insert_vertices({
        {0, 0, 0},
        {1, 0, 0},
        {2, 0, 0},
        {0, 1, 0},
        {1, 1, 0},
        {2, 1, 0},
        {0, 2, 0},
        {1, 2, 0},
        {2, 2, 0},
    });
    EXPECT_EQ(mesh->num_connected_components(), 0);

    const auto fids = mesh->insert_faces(
        {{0, 3, 1}, {1, 3, 4}, {1, 4, 2}, {2, 4, 5}, {3, 6, 4}, {4, 6, 7}, {4, 7, 5}, {5, 7, 8}});
    EXPECT_EQ(mesh->num_connected_components(), 1);

    // Insert one seam
    mesh->split_path({1, 4});
    EXPECT_EQ(mesh->num_connected_components(), 1);

    // Insert a second seam
    mesh->split_path({4, 7});
    EXPECT_EQ(mesh->num_connected_components(), 2);
}

TEST(HalfEdgeMesh, FindPath)
{
    const auto mesh = MeshType::New();
    mesh->insert_vertices({
        {0, 0, 0},
        {1, 0, 0},
        {2, 0, 0},
        {0, 1, 0},
        {1, 1, 0},
        {2, 1, 0},
        {0, 2, 0},
        {1, 2, 0},
        {2, 2, 0},
    });

    const auto fids = mesh->insert_faces(
        {{0, 3, 1}, {1, 3, 4}, {1, 4, 2}, {2, 4, 5}, {3, 6, 4}, {4, 6, 7}, {4, 7, 8}, {5, 4, 8}});

    const auto path = FindEdgePath(mesh, 1, 8);
    std::vector<std::size_t> indices;
    const std::vector<std::size_t> expected{1, 4, 8};
    indices.push_back(path[0]->vertex->idx);
    for (const auto& e : path) {
        indices.push_back(e->pair->vertex->idx);
    }
    EXPECT_EQ(indices, expected);
}

TEST(HalfEdgeMesh, FindPath_Disconnected)
{
    // Build two disconnected triangles (separate connected components)
    const auto mesh = MeshType::New();
    // CC1: vertices 0-2
    mesh->insert_vertex(0, 0, 0);
    mesh->insert_vertex(1, 0, 0);
    mesh->insert_vertex(0, 1, 0);
    mesh->insert_face(0, 2, 1);
    // CC2: vertices 3-5
    mesh->insert_vertex(5, 0, 0);
    mesh->insert_vertex(6, 0, 0);
    mesh->insert_vertex(5, 1, 0);
    mesh->insert_face(3, 5, 4);

    // Path between vertices in different CCs should be empty
    const auto path = FindEdgePath(mesh, 0, 3);
    EXPECT_TRUE(path.empty());
}

TEST(HalfEdgeMesh, IterateVertices)
{
    const auto mesh = ConstructPyramid<MeshType>();
    std::size_t count{0};
    for (const auto& v : mesh->vertices()) {
        EXPECT_EQ(v, mesh->vertex(v->idx));
        ++count;
    }
    EXPECT_EQ(count, mesh->num_vertices());
}

TEST(HalfEdgeMesh, IterateFaces)
{
    const auto mesh = ConstructPyramid<MeshType>();
    std::size_t count{0};
    for (const auto& f : mesh->faces()) {
        EXPECT_EQ(f, mesh->face(f->idx));
        ++count;
    }
    EXPECT_EQ(count, mesh->num_faces());
}

TEST(HalfEdgeMesh, IterateEdges)
{
    const auto mesh = ConstructPyramid<MeshType>();
    // Collect via mesh->edges()
    std::vector<std::size_t> edgeIdxs;
    for (const auto& e : mesh->edges()) {
        edgeIdxs.push_back(e->idxI.value());
    }
    // Collect via face-by-face iteration
    std::vector<std::size_t> expected;
    for (const auto& f : mesh->faces()) {
        for (const auto& e : *f) {
            expected.push_back(e->idxI.value());
        }
    }
    EXPECT_EQ(edgeIdxs, expected);
}

TEST(HalfEdgeMesh, IterateInteriorBoundaryPartition)
{
    const auto mesh = ConstructPyramid<MeshType>();
    std::size_t intCount{0}, bndCount{0};
    for (const auto& v : mesh->vertices_interior()) {
        EXPECT_TRUE(v->is_interior());
        ++intCount;
    }
    for (const auto& v : mesh->vertices_boundary()) {
        EXPECT_TRUE(v->is_boundary());
        ++bndCount;
    }
    EXPECT_EQ(intCount + bndCount, mesh->num_vertices());
    EXPECT_EQ(intCount, mesh->num_vertices_interior());
}

TEST(HalfEdgeMesh, IterateWheel)
{
    // Pyramid apex (vertex 3) is interior — wheel should yield all its face edges
    const auto mesh = ConstructPyramid<MeshType>();
    const auto apex = mesh->vertex(3);
    EXPECT_TRUE(apex->is_interior());

    std::vector<std::size_t> wheelIdxs;
    for (const auto& e : apex->wheel()) {
        wheelIdxs.push_back(e->idxI.value());
    }
    EXPECT_EQ(wheelIdxs.size(), 3u);  // pyramid apex touches 3 faces
    // All edges must belong to this vertex
    for (const auto& e : apex->wheel()) {
        EXPECT_EQ(e->vertex, apex);
    }
}

TEST(HalfEdgeMesh, FaceEdges)
{
    const auto mesh = ConstructPyramid<MeshType>();
    for (const auto& f : mesh->faces()) {
        // face->edges() and *face should yield the same edge sequence
        std::vector<std::size_t> viaEdges, viaStar;
        for (const auto& e : f->edges()) {
            viaEdges.push_back(e->idx);
        }
        for (const auto& e : *f) {
            viaStar.push_back(e->idx);
        }
        EXPECT_EQ(viaEdges, viaStar);
    }
}

TEST(HalfEdgeMesh, IterateEdgesCount)
{
    // edges() must yield exactly 3 * num_faces() half-edges
    const auto mesh = ConstructPyramid<MeshType>();
    std::size_t count{0};
    for (const auto& e : mesh->edges()) {
        (void)e;
        ++count;
    }
    EXPECT_EQ(count, 3 * mesh->num_faces());
}

TEST(HalfEdgeMesh, IterateWheelBoundaryVertex)
{
    // Boundary vertices still participate in faces — wheel must work and terminate
    const auto mesh = ConstructPyramid<MeshType>();
    const auto v = mesh->vertex(0);
    EXPECT_TRUE(v->is_boundary());

    std::size_t count{0};
    for (const auto& e : v->wheel()) {
        EXPECT_EQ(e->vertex, v);
        ++count;
    }
    EXPECT_EQ(count, 2u);  // vertex 0 is in exactly 2 faces in the pyramid
}

TEST(HalfEdgeMesh, IterateWheelReIterable)
{
    // Iterating a Range twice must yield identical results
    const auto mesh = ConstructPyramid<MeshType>();
    const auto apex = mesh->vertex(3);

    auto wheel = apex->wheel();
    std::vector<std::size_t> first, second;
    for (const auto& e : wheel) {
        first.push_back(e->idxI.value());
    }
    for (const auto& e : wheel) {
        second.push_back(e->idxI.value());
    }
    EXPECT_EQ(first, second);
    EXPECT_EQ(first.size(), 3u);
}

TEST(HalfEdgeMesh, IterateNoInteriorVertices)
{
    // A single triangle has no interior vertices
    const auto mesh = MeshType::New();
    mesh->insert_vertex(0, 0, 0);
    mesh->insert_vertex(1, 0, 0);
    mesh->insert_vertex(0, 1, 0);
    mesh->insert_face(0, 1, 2);

    EXPECT_TRUE(mesh->vertices_interior().empty());
    EXPECT_EQ(mesh->num_vertices_interior(), 0u);
    EXPECT_FALSE(mesh->vertices_boundary().empty());
}

TEST(HalfEdgeMesh, RangeFrontEmpty)
{
    // Range::front() and Range::empty() must behave correctly
    const auto mesh = ConstructPyramid<MeshType>();

    EXPECT_FALSE(mesh->edges().empty());
    EXPECT_NE(mesh->edges().front(), nullptr);
    EXPECT_FALSE(mesh->vertices_boundary().empty());
    EXPECT_NE(mesh->vertices_boundary().front(), nullptr);
    EXPECT_FALSE(mesh->vertices_interior().empty());
    EXPECT_NE(mesh->vertices_interior().front(), nullptr);
}

TEST(HalfEdgeMesh, Clone)
{
    const auto mesh = ConstructPyramid<MeshType>();
    const auto clone = mesh->clone();

    EXPECT_EQ(mesh->num_vertices(), clone->num_vertices());
    EXPECT_EQ(mesh->num_edges(), clone->num_edges());
    EXPECT_EQ(mesh->num_faces(), clone->num_faces());

    // Check vertices
    for (std::size_t vid = 0; vid < clone->num_vertices(); ++vid) {
        auto a = mesh->vertex(vid);
        auto b = clone->vertex(vid);
        EXPECT_NE(a, b);
        EXPECT_EQ(a->idx, b->idx);
        EXPECT_EQ(a->pos, b->pos);
        EXPECT_NE(a->edge, b->edge);
        EXPECT_EQ(a->edge->idx, b->edge->idx);
        EXPECT_NE(a->mesh, b->mesh);
    }

    // Check faces/edges
    for (std::size_t fid = 0; fid < clone->num_faces(); ++fid) {
        auto a = mesh->face(fid);
        auto b = clone->face(fid);
        EXPECT_NE(a, b);
        EXPECT_EQ(a->idx, b->idx);
        if (a->next) {
            EXPECT_NE(a->next, b->next);
            EXPECT_EQ(a->next->idx, b->next->idx);
        } else {
            EXPECT_EQ(a->next, b->next);
        }
        EXPECT_NE(a->mesh, b->mesh);

        // Check edges
        auto itA = a->begin();
        auto itB = b->begin();
        for (; itA != a->end(); ++itA, ++itB) {
            auto eA = *itA;
            auto eB = *itB;
            EXPECT_NE(eA, eB);
            EXPECT_EQ(eA->idx, eB->idx);
            EXPECT_EQ(eA->idxI, eB->idxI);

            EXPECT_NE(eA->pair, eB->pair);
            EXPECT_EQ(eA->pair->idx, eB->pair->idx);
            EXPECT_NE(eA->next, eB->next);
            EXPECT_EQ(eA->next->idx, eB->next->idx);
            EXPECT_NE(eA->prev, eB->prev);
            EXPECT_EQ(eA->prev->idx, eB->prev->idx);

            EXPECT_NE(eA->vertex, eB->vertex);
            EXPECT_EQ(eA->vertex->idx, eB->vertex->idx);

            EXPECT_NE(eA->face, eB->face);
            EXPECT_EQ(eA->face->idx, eB->face->idx);

            EXPECT_NE(eA->mesh, eB->mesh);
        }
    }
}