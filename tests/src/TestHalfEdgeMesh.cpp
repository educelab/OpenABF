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

TEST(HalfEdgeMesh, SplitMidFanAtBoundaryVertex)
{
    // 3-triangle fan around boundary vertex v0. The fan is open between v1
    // (right) and v4 (left); the existing boundary at v0 has incoming half-edge
    // 1->0 and outgoing half-edge 0->4. Splitting edge 0->2 (the middle of the
    // fan) was speculated to produce a pinched non-manifold vertex; verify
    // empirically that the partition is clean: {T_a} on one side, {T_b, T_c}
    // on the other, two connected components.
    const auto mesh = MeshType::New();
    mesh->insert_vertices({
        {0.f, 0.f, 0.f},    // v0 — fan center, on boundary
        {1.f, 0.f, 0.f},    // v1 — right fan end
        {0.5f, 1.f, 0.f},   // v2
        {-0.5f, 1.f, 0.f},  // v3
        {-1.f, 0.f, 0.f},   // v4 — left fan end
    });
    mesh->insert_faces({{0, 1, 2}, {0, 2, 3}, {0, 3, 4}});
    EXPECT_EQ(mesh->num_connected_components(), 1);

    // Split the middle edge of the fan (between T_a and T_b)
    mesh->split_path({0, 2});
    EXPECT_EQ(mesh->num_connected_components(), 2);

    // Boundary cycles should partition cleanly: one short cycle around the
    // detached T_a, one longer cycle around the T_b U T_c piece.
    const auto loops = mesh->boundaries();
    ASSERT_EQ(loops.size(), 2u);
    EXPECT_NE(loops[0].size(), loops[1].size());  // 3 vs 4

    // Every vertex on the resulting mesh must be manifold
    for (const auto& v : mesh->vertices()) {
        EXPECT_TRUE(v->is_manifold()) << "vertex " << v->idx << " is non-manifold";
    }
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

TEST(HalfEdgeMesh, ExtractConnectedComponentsSingleCCPassthrough)
{
    // Single-CC mesh: extract_connected_components returns one element whose
    // geometry and connectivity match the source and whose vertex/face maps
    // are both the identity over [0, N).
    auto mesh = MeshType::New();
    mesh->insert_vertices({{0.f, 0.f, 0.f}, {1.f, 0.f, 0.f}, {0.f, 1.f, 0.f}});
    mesh->insert_faces({{0, 1, 2}});

    auto components = mesh->extract_connected_components();
    ASSERT_EQ(components.size(), 1u);

    const auto& cc = components[0];
    EXPECT_EQ(cc.mesh->num_vertices(), 3u);
    EXPECT_EQ(cc.mesh->num_faces(), 1u);
    EXPECT_EQ(cc.vertex_map, (std::vector<std::size_t>{0, 1, 2}));
    EXPECT_EQ(cc.face_map, (std::vector<std::size_t>{0}));

    for (std::size_t i = 0; i < cc.mesh->num_vertices(); ++i) {
        EXPECT_EQ(cc.mesh->vertex(i)->pos, mesh->vertex(cc.vertex_map[i])->pos);
    }
}

TEST(HalfEdgeMesh, ExtractConnectedComponentsTwoCCs)
{
    // Two disjoint triangles -> two extracted meshes with disjoint vertex
    // maps covering every original vertex exactly once, and face maps
    // pointing back to the source face indices.
    auto mesh = MeshType::New();
    mesh->insert_vertices({
        {0.f, 0.f, 0.f},
        {1.f, 0.f, 0.f},
        {0.f, 1.f, 0.f},
        {5.f, 5.f, 0.f},
        {6.f, 5.f, 0.f},
        {5.f, 6.f, 0.f},
    });
    mesh->insert_faces({{0, 1, 2}, {3, 4, 5}});
    EXPECT_EQ(mesh->num_connected_components(), 2u);

    auto components = mesh->extract_connected_components();
    ASSERT_EQ(components.size(), 2u);

    std::vector<std::size_t> originalsCovered;
    std::vector<std::size_t> facesCovered;
    for (const auto& cc : components) {
        EXPECT_EQ(cc.mesh->num_vertices(), 3u);
        EXPECT_EQ(cc.mesh->num_faces(), 1u);
        ASSERT_EQ(cc.face_map.size(), 1u);
        facesCovered.push_back(cc.face_map[0]);

        for (std::size_t i = 0; i < cc.mesh->num_vertices(); ++i) {
            EXPECT_EQ(cc.mesh->vertex(i)->idx, i);
            EXPECT_EQ(cc.mesh->vertex(i)->pos, mesh->vertex(cc.vertex_map[i])->pos);
            originalsCovered.push_back(cc.vertex_map[i]);
        }
    }

    std::sort(originalsCovered.begin(), originalsCovered.end());
    EXPECT_EQ(originalsCovered, (std::vector<std::size_t>{0, 1, 2, 3, 4, 5}));
    std::sort(facesCovered.begin(), facesCovered.end());
    EXPECT_EQ(facesCovered, (std::vector<std::size_t>{0, 1}));
}

namespace
{
// Custom edge traits with a field that has no geometric meaning, so the only
// way it can survive extraction is if clone_face_ explicitly copied it.
// Inherits DefaultEdgeTraits because insert_face_ runs ComputeFaceAngles
// which needs `alpha`.
struct CustomEdgeTraits : traits::DefaultEdgeTraits<float> {
    int label{0};
};
using LabeledMesh = HalfEdgeMesh<float, 3, traits::DefaultVertexTraits<float>, CustomEdgeTraits,
                                 traits::DefaultFaceTraits<float>>;
}  // namespace

TEST(HalfEdgeMesh, ExtractPreservesCustomEdgeTraits)
{
    // Build a two-CC mesh with a custom EdgeTraits struct. Set a unique label
    // on every half-edge before extraction; verify extracted edges carry the
    // same labels. This is the actual test for the "deep copy preserves
    // traits" guarantee — geometric alpha is recomputed by insert_face_ and
    // would survive even without explicit copying.
    auto mesh = LabeledMesh::New();
    mesh->insert_vertices({
        {0.f, 0.f, 0.f},
        {1.f, 0.f, 0.f},
        {0.f, 1.f, 0.f},
        {5.f, 5.f, 0.f},
        {6.f, 5.f, 0.f},
        {5.f, 6.f, 0.f},
    });
    mesh->insert_faces({{0, 1, 2}, {3, 4, 5}});

    // Assign a label per (face_idx, edge_in_face_idx) pair so we can match
    // them up after extraction.
    int next = 1;
    std::vector<std::vector<int>> srcLabels;
    for (const auto& f : mesh->faces()) {
        std::vector<int> faceLabels;
        for (const auto& e : *f) {
            e->label = next;
            faceLabels.push_back(next);
            ++next;
        }
        srcLabels.push_back(std::move(faceLabels));
    }

    auto components = mesh->extract_connected_components();
    ASSERT_EQ(components.size(), 2u);

    for (std::size_t c = 0; c < components.size(); ++c) {
        auto& sub = components[c].mesh;
        ASSERT_EQ(sub->num_faces(), 1u);
        std::vector<int> subLabels;
        for (const auto& e : *sub->face(0)) {
            subLabels.push_back(e->label);
        }
        EXPECT_EQ(subLabels, srcLabels[c])
            << "edge labels not preserved on extracted component " << c;
    }
}

TEST(HalfEdgeMesh, ExtractFaceMapAfterSplit)
{
    // Realistic case: 3x3 grid torn down the middle. Extracted face_maps
    // partition the source face indices [0, num_faces) without overlap.
    auto mesh = MeshType::New();
    mesh->insert_vertices({
        {0.f, 0.f, 0.f},
        {1.f, 0.f, 0.f},
        {2.f, 0.f, 0.f},
        {0.f, 1.f, 0.f},
        {1.f, 1.f, 0.f},
        {2.f, 1.f, 0.f},
        {0.f, 2.f, 0.f},
        {1.f, 2.f, 0.f},
        {2.f, 2.f, 0.f},
    });
    mesh->insert_faces({
        {0, 3, 1},
        {1, 3, 4},
        {1, 4, 2},
        {2, 4, 5},
        {3, 6, 4},
        {4, 6, 7},
        {4, 7, 5},
        {5, 7, 8},
    });
    mesh->split_path({1, 4});
    mesh->split_path({4, 7});
    ASSERT_EQ(mesh->num_connected_components(), 2u);

    auto components = mesh->extract_connected_components();
    ASSERT_EQ(components.size(), 2u);

    std::vector<std::size_t> allFaces;
    for (const auto& cc : components) {
        EXPECT_EQ(cc.face_map.size(), cc.mesh->num_faces());
        for (std::size_t i = 0; i < cc.mesh->num_faces(); ++i) {
            // Vertex positions on the extracted face match the source face
            // they map back to.
            auto srcFace = mesh->face(cc.face_map[i]);
            auto srcIt = srcFace->begin();
            auto subIt = cc.mesh->face(i)->begin();
            while (srcIt != srcFace->end()) {
                EXPECT_EQ((*srcIt)->vertex->pos, (*subIt)->vertex->pos);
                ++srcIt;
                ++subIt;
            }
            allFaces.push_back(cc.face_map[i]);
        }
    }

    std::sort(allFaces.begin(), allFaces.end());
    EXPECT_EQ(allFaces, (std::vector<std::size_t>{0, 1, 2, 3, 4, 5, 6, 7}));
}