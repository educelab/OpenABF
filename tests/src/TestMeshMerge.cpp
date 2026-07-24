#include <set>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "OpenABF/OpenABF.hpp"

using Mesh = OpenABF::HalfEdgeMesh<float>;
using OpenABF::MergedMesh;
using OpenABF::MergeMeshes;

namespace
{

/** @brief Build a w x h rectangle chart (two triangles) translated to (ox,oy) */
auto MakeRectChart(float w, float h, float ox, float oy) -> Mesh::Pointer
{
    auto m = Mesh::New();
    m->insert_vertices(
        {{ox, oy, 0.f}, {ox + w, oy, 0.f}, {ox + w, oy + h, 0.f}, {ox, oy + h, 0.f}});
    m->insert_faces({{0, 1, 2}, {0, 2, 3}});
    return m;
}

}  // namespace

TEST(MeshMerge, EmptyListGivesEmptyMesh)
{
    std::vector<Mesh::Pointer> charts;
    auto merged = MergeMeshes<Mesh>(charts);
    ASSERT_NE(merged.mesh, nullptr);
    EXPECT_EQ(merged.mesh->num_vertices(), 0u);
    EXPECT_EQ(merged.mesh->num_faces(), 0u);
    EXPECT_TRUE(merged.vertex_source.empty());
    EXPECT_TRUE(merged.face_source.empty());
}

TEST(MeshMerge, NullChartThrows)
{
    std::vector<Mesh::Pointer> charts{MakeRectChart(1.f, 1.f, 0.f, 0.f), nullptr};
    EXPECT_THROW(MergeMeshes<Mesh>(charts), std::invalid_argument);
}

TEST(MeshMerge, EmptyChartThrows)
{
    std::vector<Mesh::Pointer> charts{MakeRectChart(1.f, 1.f, 0.f, 0.f), Mesh::New()};
    EXPECT_THROW(MergeMeshes<Mesh>(charts), std::invalid_argument);
}

TEST(MeshMerge, ConcatenatesVerticesAndFaces)
{
    std::vector<Mesh::Pointer> charts{MakeRectChart(1.f, 1.f, 0.f, 0.f),
                                      MakeRectChart(1.f, 1.f, 5.f, 0.f)};
    // MeshType is deduced from the mesh vector here; the explicit spelling used
    // elsewhere in this file must keep working too.
    auto merged = MergeMeshes(charts);
    EXPECT_EQ(merged.mesh->num_vertices(), 8u);
    EXPECT_EQ(merged.mesh->num_faces(), 4u);
    EXPECT_EQ(merged.vertex_source.size(), 8u);
    EXPECT_EQ(merged.face_source.size(), 4u);
}

TEST(MeshMerge, VertexSourcePreservesPositionsAndProvenance)
{
    std::vector<Mesh::Pointer> charts{MakeRectChart(1.f, 1.f, 0.f, 0.f),
                                      MakeRectChart(2.f, 3.f, 5.f, 0.f)};
    auto merged = MergeMeshes<Mesh>(charts);

    ASSERT_EQ(merged.vertex_source.size(), merged.mesh->num_vertices());
    for (std::size_t i = 0; i < merged.mesh->num_vertices(); ++i) {
        auto [chart, sub] = merged.vertex_source[i];
        ASSERT_LT(chart, charts.size());
        ASSERT_LT(sub, charts[chart]->num_vertices());
        const auto& mergedV = merged.mesh->vertex(i);
        const auto& srcV = charts[chart]->vertex(sub);
        EXPECT_FLOAT_EQ(mergedV->pos[0], srcV->pos[0]) << "vertex " << i;
        EXPECT_FLOAT_EQ(mergedV->pos[1], srcV->pos[1]) << "vertex " << i;
    }
}

TEST(MeshMerge, FaceSourceProvenanceIsValid)
{
    std::vector<Mesh::Pointer> charts{MakeRectChart(1.f, 1.f, 0.f, 0.f),
                                      MakeRectChart(1.f, 1.f, 5.f, 0.f)};
    auto merged = MergeMeshes<Mesh>(charts);

    ASSERT_EQ(merged.face_source.size(), merged.mesh->num_faces());
    // Each (chart, sub face) provenance is unique and in range.
    std::set<std::pair<std::size_t, std::size_t>> seen;
    for (const auto& [chart, sub] : merged.face_source) {
        ASSERT_LT(chart, charts.size());
        ASSERT_LT(sub, charts[chart]->num_faces());
        EXPECT_TRUE(seen.emplace(chart, sub).second) << "duplicate face provenance";
    }
}

// Uniqueness is not enough: face_source[i] must describe merged face i. Charts
// are placed far apart so each source face has a distinct centroid.
TEST(MeshMerge, FaceSourceIndexAlignsWithMergedFace)
{
    std::vector<Mesh::Pointer> charts{MakeRectChart(1.f, 1.f, 0.f, 0.f),
                                      MakeRectChart(2.f, 3.f, 10.f, 20.f)};
    auto merged = MergeMeshes<Mesh>(charts);
    ASSERT_EQ(merged.face_source.size(), merged.mesh->num_faces());

    auto centroid = [](const auto& face) {
        float cx{0.f};
        float cy{0.f};
        float n{0.f};
        for (const auto& edge : *face) {
            cx += edge->vertex->pos[0];
            cy += edge->vertex->pos[1];
            n += 1.f;
        }
        return std::pair<float, float>{cx / n, cy / n};
    };

    for (std::size_t i = 0; i < merged.mesh->num_faces(); ++i) {
        const auto [chart, sub] = merged.face_source[i];
        ASSERT_LT(chart, charts.size());
        ASSERT_LT(sub, charts[chart]->num_faces());
        const auto [mx, my] = centroid(merged.mesh->face(i));
        const auto [sx, sy] = centroid(charts[chart]->face(sub));
        EXPECT_FLOAT_EQ(mx, sx) << "merged face " << i;
        EXPECT_FLOAT_EQ(my, sy) << "merged face " << i;
    }
}

// Merged vertices must reference half-edges of the merged mesh, never of a
// source mesh: Vertex::is_boundary(), wheel(), and is_unreferenced() all read
// `edge` directly, so a carried-over pointer would silently report the source
// chart's neighbourhood.
TEST(MeshMerge, MergedVerticesReferenceMergedMeshEdges)
{
    std::vector<Mesh::Pointer> charts{MakeRectChart(1.f, 1.f, 0.f, 0.f),
                                      MakeRectChart(1.f, 1.f, 5.f, 0.f)};
    auto merged = MergeMeshes<Mesh>(charts);
    for (const auto& v : merged.mesh->vertices()) {
        ASSERT_NE(v->edge, nullptr) << "vertex " << v->idx << " has no edge";
        EXPECT_EQ(v->edge->mesh, merged.mesh.get())
            << "vertex " << v->idx << " references an edge outside the merged mesh";
        EXPECT_EQ(v->mesh, merged.mesh.get()) << "vertex " << v->idx;
    }
    // Every vertex of these two disjoint quads is on a boundary; this traverses
    // via vertex->edge, so it only holds if the pointers are rebound correctly.
    for (const auto& v : merged.mesh->vertices()) {
        EXPECT_TRUE(v->is_boundary()) << "vertex " << v->idx;
    }
}

// Composing the merge maps with extract's back-maps must recover original
// (torn-mesh) identity: this is the full chain merged -> chart -> M'.
TEST(MeshMerge, RoundTripExtractMergeRecoversOriginalIdentity)
{
    using ABF = OpenABF::ABFPlusPlus<float>;
    using LSCM = OpenABF::AngleBasedLSCM<float, ABF::Mesh>;

    auto mesh = ABF::Mesh::New();
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
    mesh->split_path({1, 4, 7});
    auto ccs = mesh->extract_connected_components();
    ASSERT_EQ(ccs.size(), 2u);

    std::vector<ABF::Mesh::Pointer> chartMeshes;
    std::size_t totalVerts = 0;
    std::size_t totalFaces = 0;
    for (auto& cc : ccs) {
        LSCM::Compute(cc.mesh);
        chartMeshes.push_back(cc.mesh);
        totalVerts += cc.mesh->num_vertices();
        totalFaces += cc.mesh->num_faces();
    }

    auto merged = MergeMeshes<ABF::Mesh>(chartMeshes);
    EXPECT_EQ(merged.mesh->num_vertices(), totalVerts);
    EXPECT_EQ(merged.mesh->num_faces(), totalFaces);

    // merged face -> (chart, sub face) -> M' face via extract face_map.
    // Every original (torn-mesh) face must be hit exactly once.
    std::set<std::size_t> origFaces;
    for (const auto& [chart, sub] : merged.face_source) {
        auto origFace = ccs[chart].face_map[sub];
        EXPECT_TRUE(origFaces.insert(origFace).second) << "duplicate original face " << origFace;
    }
    EXPECT_EQ(origFaces.size(), 8u);  // all original faces recovered
}
