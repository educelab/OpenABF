#include <gtest/gtest.h>

#include "OpenABF/OpenABF.hpp"

using namespace OpenABF;

// Use the ABF-traited mesh so the same MeshType works for both extraction
// tests and the end-to-end ParameterizeConnectedComponents test (which needs
// the ABF::Mesh edge/face traits for ABFPlusPlus + AngleBasedLSCM).
using Optimizer = ABFPlusPlus<float>;
using MeshType = Optimizer::Mesh;
using Param = AngleBasedLSCM<float, MeshType>;

namespace
{
auto MakeTwoCCMesh() -> MeshType::Pointer
{
    // Two disjoint triangles, no shared vertices.
    auto mesh = MeshType::New();
    mesh->insert_vertices({
        // CC0 — original-idx 0..2
        {0.f, 0.f, 0.f},
        {1.f, 0.f, 0.f},
        {0.f, 1.f, 0.f},
        // CC1 — original-idx 3..5
        {5.f, 5.f, 0.f},
        {6.f, 5.f, 0.f},
        {5.f, 6.f, 0.f},
    });
    mesh->insert_faces({{0, 1, 2}, {3, 4, 5}});
    return mesh;
}
}  // namespace

TEST(HalfEdgeMeshUtils, ExtractSingleCCPassthrough)
{
    // Single-CC mesh: ExtractConnectedComponents returns one element whose
    // geometry and connectivity match the source.
    auto mesh = MeshType::New();
    mesh->insert_vertices({{0.f, 0.f, 0.f}, {1.f, 0.f, 0.f}, {0.f, 1.f, 0.f}});
    mesh->insert_faces({{0, 1, 2}});

    auto components = ExtractConnectedComponents<MeshType>(mesh);
    ASSERT_EQ(components.size(), 1u);

    const auto& [sub, backMap] = components[0];
    EXPECT_EQ(sub->num_vertices(), 3u);
    EXPECT_EQ(sub->num_faces(), 1u);
    EXPECT_EQ(backMap.size(), 3u);

    // Each extracted vertex's pos matches the original at backMap[i]
    for (std::size_t i = 0; i < sub->num_vertices(); ++i) {
        EXPECT_EQ(sub->vertex(i)->pos, mesh->vertex(backMap[i])->pos);
    }
}

TEST(HalfEdgeMeshUtils, ExtractTwoCCs)
{
    // Two disjoint triangles -> two extracted meshes with disjoint back-maps
    auto mesh = MakeTwoCCMesh();
    EXPECT_EQ(mesh->num_connected_components(), 2u);

    auto components = ExtractConnectedComponents<MeshType>(mesh);
    ASSERT_EQ(components.size(), 2u);

    std::vector<std::size_t> originalsCovered;
    for (const auto& [sub, backMap] : components) {
        EXPECT_EQ(sub->num_vertices(), 3u);
        EXPECT_EQ(sub->num_faces(), 1u);

        // Indices in extracted mesh are 0..N-1
        for (std::size_t i = 0; i < sub->num_vertices(); ++i) {
            EXPECT_EQ(sub->vertex(i)->idx, i);
            EXPECT_EQ(sub->vertex(i)->pos, mesh->vertex(backMap[i])->pos);
            originalsCovered.push_back(backMap[i]);
        }
    }

    // Every original vertex appears in exactly one extracted mesh
    std::sort(originalsCovered.begin(), originalsCovered.end());
    EXPECT_EQ(originalsCovered, (std::vector<std::size_t>{0, 1, 2, 3, 4, 5}));
}

TEST(HalfEdgeMeshUtils, ExtractPreservesEdgeAlphaTraits)
{
    // Set face angle traits on the source, extract, verify traits are copied
    // onto the extracted mesh's edges.
    auto mesh = MakeTwoCCMesh();
    ComputeMeshAngles(mesh);

    // Snapshot the (face_idx, edge-in-face, alpha) tuples for later comparison
    std::vector<std::vector<float>> srcAlphas;
    for (const auto& f : mesh->faces()) {
        std::vector<float> faceAlphas;
        for (const auto& e : *f) {
            faceAlphas.push_back(e->alpha);
            EXPECT_GT(e->alpha, 0.f);  // sanity: angles really did compute
        }
        srcAlphas.push_back(std::move(faceAlphas));
    }

    auto components = ExtractConnectedComponents<MeshType>(mesh);
    ASSERT_EQ(components.size(), 2u);

    // Extracted faces line up with components in insertion order (CC discovery
    // follows face insertion order in num_connected_components()).
    for (std::size_t c = 0; c < components.size(); ++c) {
        auto& sub = components[c].first;
        ASSERT_EQ(sub->num_faces(), 1u);
        std::vector<float> subAlphas;
        for (const auto& e : *sub->face(0)) {
            subAlphas.push_back(e->alpha);
        }
        EXPECT_EQ(subAlphas, srcAlphas[c]);
    }
}

TEST(HalfEdgeMeshUtils, ParameterizeWritesUVBackToOriginal)
{
    // Two-CC mesh with 3D positions; ParameterizeConnectedComponents should
    // flatten each component and write UV coords back into the source's
    // vertex positions.
    auto mesh = MakeTwoCCMesh();

    // Snapshot original 3D positions so we can detect that they were rewritten
    std::vector<Vec<float, 3>> before;
    for (const auto& v : mesh->vertices()) {
        before.push_back(v->pos);
    }

    EXPECT_NO_THROW((ParameterizeConnectedComponents<Optimizer, Param>(mesh)));

    // At least one vertex per component should now have its z coordinate
    // zeroed (LSCM lives in 2D and parameterizers in this library typically
    // map onto the XY plane).
    bool cc0Changed = false;
    bool cc1Changed = false;
    for (std::size_t i = 0; i < mesh->num_vertices(); ++i) {
        if (mesh->vertex(i)->pos != before[i]) {
            if (i < 3) {
                cc0Changed = true;
            } else {
                cc1Changed = true;
            }
        }
    }
    EXPECT_TRUE(cc0Changed);
    EXPECT_TRUE(cc1Changed);
}
