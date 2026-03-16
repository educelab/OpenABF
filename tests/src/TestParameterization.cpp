#include <chrono>
#include <cmath>
#include <iostream>
#include <gtest/gtest.h>

#include "OpenABF/OpenABF.hpp"
#include "Utils.hpp"

using namespace OpenABF;
using namespace OpenABF::tests;

TEST(Parameterization, AngledBasedLSCM)
{
    using LSCM = AngleBasedLSCM<float>;

    auto mesh = ConstructPyramid<LSCM::Mesh>();
    LSCM::Compute(mesh);

    const std::vector expected{Vec3f{0, 0, 0}, Vec3f{2, 0, 0}, Vec3f{1, 1.0392305, 0},
                               Vec3f{1, 0.34641013, 0}};
    for (auto v = 0; v < mesh->num_vertices(); ++v) {
        const auto& vv = mesh->vertex(v);
        const auto& ve = expected[v];
        for (auto i = 0; i < 3; i++) {
            EXPECT_FLOAT_EQ(vv->pos[i], ve[i]);
        }
    }
}

TEST(Parameterizations, ABF)
{
    using ABF = ABF<float>;
    using LSCM = AngleBasedLSCM<float, ABF::Mesh>;

    auto mesh = ConstructPyramid<ABF::Mesh>();
    ABF::Compute(mesh);
    LSCM::Compute(mesh);

    const std::vector expected{Vec3f{0, 0, 0}, Vec3f{2, 0, 0}, Vec3f{1, 1.7320509, 0},
                               Vec3f{1, 0.5773503, 0}};
    for (auto v = 0; v < mesh->num_vertices(); ++v) {
        const auto& vv = mesh->vertex(v);
        const auto& ve = expected[v];
        for (auto i = 0; i < 3; i++) {
            EXPECT_FLOAT_EQ(vv->pos[i], ve[i]);
        }
    }
}

TEST(Parameterizations, ABFPlusPlus)
{
    using ABF = ABFPlusPlus<float>;
    using LSCM = AngleBasedLSCM<float, ABF::Mesh>;

    auto mesh = ConstructPyramid<ABF::Mesh>();
    ABF::Compute(mesh);
    LSCM::Compute(mesh);

    const std::vector expected{Vec3f{0, 0, 0}, Vec3f{2, 0, 0}, Vec3f{1, 1.7320509, 0},
                               Vec3f{1, 0.5773503, 0}};
    for (auto v = 0; v < mesh->num_vertices(); ++v) {
        const auto& vv = mesh->vertex(v);
        const auto& ve = expected[v];
        for (auto i = 0; i < 3; i++) {
            EXPECT_FLOAT_EQ(vv->pos[i], ve[i]);
        }
    }
}

TEST(Parameterizations, ABFPlusPlus_Double)
{
    // Smoke test confirming ABFPlusPlus<double> produces the same result as
    // ABFPlusPlus<float> to within double precision (catches the 1.F / T(1) bug)
    using ABF = ABFPlusPlus<double>;
    using LSCM = AngleBasedLSCM<double, ABF::Mesh>;

    auto mesh = ConstructPyramid<ABF::Mesh>();
    ABF::Compute(mesh);
    LSCM::Compute(mesh);

    using Vec3d = Vec<double, 3>;
    const std::vector expected{Vec3d{0, 0, 0}, Vec3d{2, 0, 0}, Vec3d{1, 1.7320508075688772, 0},
                               Vec3d{1, 0.5773502691896258, 0}};
    for (auto v = 0; v < mesh->num_vertices(); ++v) {
        const auto& vv = mesh->vertex(v);
        const auto& ve = expected[v];
        for (auto i = 0; i < 3; i++) {
            EXPECT_DOUBLE_EQ(vv->pos[i], ve[i]);
        }
    }
}

TEST(Parameterizations, ABF_Double)
{
    using ABFType = ABF<double>;
    using LSCM = AngleBasedLSCM<double, ABFType::Mesh>;

    auto mesh = ConstructPyramid<ABFType::Mesh>();
    ABFType::Compute(mesh);
    LSCM::Compute(mesh);

    using Vec3d = Vec<double, 3>;
    const std::vector expected{Vec3d{0, 0, 0}, Vec3d{2, 0, 0}, Vec3d{1, 1.7320508075688772, 0},
                               Vec3d{1, 0.5773502691896258, 0}};
    for (auto v = 0; v < mesh->num_vertices(); ++v) {
        const auto& vv = mesh->vertex(v);
        const auto& ve = expected[v];
        for (auto i = 0; i < 3; i++) {
            EXPECT_DOUBLE_EQ(vv->pos[i], ve[i]);
        }
    }
}

TEST(Parameterizations, AngleBasedLSCM_Double)
{
    using LSCM = AngleBasedLSCM<double>;

    auto mesh = ConstructPyramid<LSCM::Mesh>();
    LSCM::Compute(mesh);

    using Vec3d = Vec<double, 3>;
    const std::vector expected{Vec3d{0, 0, 0}, Vec3d{2, 0, 0},
                               Vec3d{1.0000000000000002, 1.0392304845413267, 0},
                               Vec3d{1, 0.3464101615137756, 0}};
    for (auto v = 0; v < mesh->num_vertices(); ++v) {
        const auto& vv = mesh->vertex(v);
        const auto& ve = expected[v];
        for (auto i = 0; i < 3; i++) {
            EXPECT_DOUBLE_EQ(vv->pos[i], ve[i]);
        }
    }
}

TEST(Parameterizations, ABF_MultiInterior)
{
    // 4×4 vertex grid → 4 interior vertices; exercises the multi-interior code path
    using ABFType = ABF<float>;
    using LSCM = AngleBasedLSCM<float, ABFType::Mesh>;

    auto mesh = ConstructGrid<ABFType::Mesh>(4, 4);
    ASSERT_EQ(mesh->num_vertices_interior(), 4u);

    ABFType::Compute(mesh);
    LSCM::Compute(mesh);

    // After LSCM, all UV z-coordinates must be zero and x/y must be finite
    for (std::size_t v = 0; v < mesh->num_vertices(); ++v) {
        const auto& pos = mesh->vertex(v)->pos;
        EXPECT_TRUE(std::isfinite(pos[0])) << "vertex " << v << " x is not finite";
        EXPECT_TRUE(std::isfinite(pos[1])) << "vertex " << v << " y is not finite";
        EXPECT_FLOAT_EQ(pos[2], 0.f);
    }
}

TEST(Parameterizations, ABFPlusPlus_MultiInterior)
{
    // 4×4 vertex grid → 4 interior vertices
    using ABFType = ABFPlusPlus<float>;
    using LSCM = AngleBasedLSCM<float, ABFType::Mesh>;

    auto mesh = ConstructGrid<ABFType::Mesh>(4, 4);
    ASSERT_EQ(mesh->num_vertices_interior(), 4u);

    ABFType::Compute(mesh);
    LSCM::Compute(mesh);

    for (std::size_t v = 0; v < mesh->num_vertices(); ++v) {
        const auto& pos = mesh->vertex(v)->pos;
        EXPECT_TRUE(std::isfinite(pos[0])) << "vertex " << v << " x is not finite";
        EXPECT_TRUE(std::isfinite(pos[1])) << "vertex " << v << " y is not finite";
        EXPECT_FLOAT_EQ(pos[2], 0.f);
    }
}

TEST(Parameterizations, ABF_NoInteriorVertices)
{
    // Single triangle: 0 interior vertices — solver must not crash
    using ABFType = ABF<float>;
    using LSCM = AngleBasedLSCM<float, ABFType::Mesh>;

    auto mesh = ABFType::Mesh::New();
    mesh->insert_vertex(0, 0, 0);
    mesh->insert_vertex(1, 0, 0);
    mesh->insert_vertex(0, 1, 0);
    mesh->insert_face(0, 2, 1);

    ASSERT_EQ(mesh->num_vertices_interior(), 0u);
    EXPECT_NO_THROW(ABFType::Compute(mesh));
    EXPECT_NO_THROW(LSCM::Compute(mesh));
}

TEST(Parameterizations, ABF_MaxIters)
{
    // Instance API: limit to 1 iteration and confirm the solver respects it
    using ABFType = ABF<float>;

    auto mesh = ConstructPyramid<ABFType::Mesh>();
    ABFType abf;
    abf.setMaxIterations(1);
    abf.compute(mesh);

    EXPECT_EQ(abf.iterations(), 1u);
}

TEST(Parameterizations, ABFPlusPlus_TightThreshold)
{
    // Tighter threshold should produce a lower final gradient than the default
    using ABFType = ABFPlusPlus<float>;

    auto mesh_default = ConstructPyramid<ABFType::Mesh>();
    ABFType abf_default;
    abf_default.compute(mesh_default);

    auto mesh_tight = ConstructPyramid<ABFType::Mesh>();
    ABFType abf_tight;
    abf_tight.setGradientThreshold(1e-6f);
    abf_tight.setMaxIterations(100);
    abf_tight.compute(mesh_tight);

    EXPECT_LE(abf_tight.gradient(), abf_default.gradient());
}

TEST(Parameterizations, ABFPlusPlus_LooseThreshold)
{
    // A threshold larger than the initial gradient causes zero solver iterations
    using ABFType = ABFPlusPlus<float>;

    auto mesh = ConstructPyramid<ABFType::Mesh>();
    ABFType abf;
    abf.setGradientThreshold(1e6f);
    abf.compute(mesh);

    EXPECT_EQ(abf.iterations(), 0u);
}

TEST(Parameterization, AngleBasedLSCM_ExplicitPins)
{
    // Explicit pins matching the default selection must produce identical results
    using LSCM = AngleBasedLSCM<float>;

    auto mesh = ConstructPyramid<LSCM::Mesh>();
    LSCM::Compute(mesh, 0, 1);

    const std::vector expected{Vec3f{0, 0, 0}, Vec3f{2, 0, 0}, Vec3f{1, 1.0392305, 0},
                               Vec3f{1, 0.34641013, 0}};
    for (auto v = 0; v < mesh->num_vertices(); ++v) {
        const auto& vv = mesh->vertex(v);
        const auto& ve = expected[v];
        for (auto i = 0; i < 3; i++) {
            EXPECT_FLOAT_EQ(vv->pos[i], ve[i]);
        }
    }
}

TEST(Parameterization, AngleBasedLSCM_ExplicitPins_Reversed)
{
    // Swapping the pin pair must change which vertex lands at the origin
    using LSCM = AngleBasedLSCM<float>;

    auto mesh = ConstructPyramid<LSCM::Mesh>();
    LSCM::Compute(mesh, 1, 0);

    // p0 = vertex 1 → placed at {0, 0, 0}
    // p1 = vertex 0 → edge from v1 to v0 is {-1,0,0}; max axis is Y (value 0),
    //   copysign(dist=2, 0) = +2  →  p1 placed at {0, 2, 0}
    EXPECT_FLOAT_EQ(mesh->vertex(1)->pos[0], 0.f);
    EXPECT_FLOAT_EQ(mesh->vertex(1)->pos[1], 0.f);
    EXPECT_FLOAT_EQ(mesh->vertex(0)->pos[0], 0.f);
    EXPECT_FLOAT_EQ(mesh->vertex(0)->pos[1], 2.f);
}

TEST(Parameterization, AngleBasedLSCM_SetPinnedVertices)
{
    // Instance API: setPinnedVertices selects the same pins as the static overload
    using LSCM = AngleBasedLSCM<float>;

    auto mesh_static = ConstructPyramid<LSCM::Mesh>();
    LSCM::Compute(mesh_static, 0, 1);

    auto mesh_instance = ConstructPyramid<LSCM::Mesh>();
    LSCM lscm;
    lscm.setPinnedVertices(0, 1);
    lscm.compute(mesh_instance);

    for (auto v = 0; v < mesh_static->num_vertices(); ++v) {
        const auto& vs = mesh_static->vertex(v)->pos;
        const auto& vi = mesh_instance->vertex(v)->pos;
        for (auto i = 0; i < 3; i++) {
            EXPECT_FLOAT_EQ(vi[i], vs[i]);
        }
    }
}

TEST(Fixtures, Hemisphere)
{
    using M = HalfEdgeMesh<float>;
    auto mesh = ConstructHemisphere<M>(12, 24);
    EXPECT_EQ(mesh->num_vertices(), 289u);
    EXPECT_EQ(mesh->num_faces(), 552u);
    EXPECT_GT(mesh->num_vertices_interior(), 0u);
}

TEST(Fixtures, WavySurface)
{
    using M = HalfEdgeMesh<float>;
    auto mesh = ConstructWavySurface<M>(20, 20);
    EXPECT_EQ(mesh->num_vertices(), 400u);
    EXPECT_EQ(mesh->num_faces(), 722u);
    EXPECT_GT(mesh->num_vertices_interior(), 0u);
}

// ============================================================================
// HLSCM tests (F1)
// ============================================================================

TEST(HLSCM, Pyramid)
{
    // Single-level fallback: HLSCM must produce identical output to LSCM
    using HLSCM = HierarchicalLSCM<float>;
    using LSCM = AngleBasedLSCM<float>;

    auto mesh_lscm = ConstructPyramid<LSCM::Mesh>();
    LSCM::Compute(mesh_lscm);

    auto mesh_hlscm = ConstructPyramid<HLSCM::Mesh>();
    HLSCM::Compute(mesh_hlscm);

    for (std::size_t v = 0; v < mesh_lscm->num_vertices(); ++v) {
        for (auto i = 0; i < 3; i++) {
            EXPECT_FLOAT_EQ(mesh_hlscm->vertex(v)->pos[i], mesh_lscm->vertex(v)->pos[i])
                << "vertex " << v << " component " << i;
        }
    }
}

TEST(HLSCM, ExplicitPins)
{
    using HLSCM = HierarchicalLSCM<float>;

    auto mesh_auto = ConstructPyramid<HLSCM::Mesh>();
    HLSCM::Compute(mesh_auto);

    auto mesh_explicit = ConstructPyramid<HLSCM::Mesh>();
    HLSCM::Compute(mesh_explicit, 0, 1);

    for (std::size_t v = 0; v < mesh_auto->num_vertices(); ++v) {
        for (auto i = 0; i < 3; i++) {
            EXPECT_FLOAT_EQ(mesh_explicit->vertex(v)->pos[i], mesh_auto->vertex(v)->pos[i]);
        }
    }
}

TEST(HLSCM, SetPinnedVertices)
{
    using HLSCM = HierarchicalLSCM<float>;

    auto mesh_static = ConstructPyramid<HLSCM::Mesh>();
    HLSCM::Compute(mesh_static, 0, 1);

    auto mesh_instance = ConstructPyramid<HLSCM::Mesh>();
    HLSCM hlscm;
    hlscm.setPinnedVertices(0, 1);
    hlscm.compute(mesh_instance);

    for (std::size_t v = 0; v < mesh_static->num_vertices(); ++v) {
        for (auto i = 0; i < 3; i++) {
            EXPECT_FLOAT_EQ(mesh_instance->vertex(v)->pos[i], mesh_static->vertex(v)->pos[i]);
        }
    }
}

TEST(HLSCM, Double)
{
    using HLSCM = HierarchicalLSCM<double>;
    using LSCM = AngleBasedLSCM<double>;

    auto mesh_lscm = ConstructPyramid<LSCM::Mesh>();
    LSCM::Compute(mesh_lscm);

    auto mesh_hlscm = ConstructPyramid<HLSCM::Mesh>();
    HLSCM::Compute(mesh_hlscm);

    for (std::size_t v = 0; v < mesh_lscm->num_vertices(); ++v) {
        for (auto i = 0; i < 3; i++) {
            EXPECT_DOUBLE_EQ(mesh_hlscm->vertex(v)->pos[i], mesh_lscm->vertex(v)->pos[i]);
        }
    }
}

TEST(HLSCM, ABFPlusPlus)
{
    using ABFType = ABFPlusPlus<float>;
    using HLSCM = HierarchicalLSCM<float, ABFType::Mesh>;

    auto mesh = ConstructPyramid<ABFType::Mesh>();
    ABFType::Compute(mesh);
    HLSCM::Compute(mesh);

    for (std::size_t v = 0; v < mesh->num_vertices(); ++v) {
        const auto& pos = mesh->vertex(v)->pos;
        EXPECT_TRUE(std::isfinite(pos[0])) << "vertex " << v;
        EXPECT_TRUE(std::isfinite(pos[1])) << "vertex " << v;
        EXPECT_FLOAT_EQ(pos[2], 0.f);
    }
}

TEST(HLSCM, Hemisphere)
{
    // Non-trivial curvature: verify valid parameterization and no flipped triangles
    using HLSCM = HierarchicalLSCM<float>;
    auto mesh = ConstructHemisphere<HLSCM::Mesh>(12, 24);
    HLSCM::Compute(mesh);

    for (std::size_t v = 0; v < mesh->num_vertices(); ++v) {
        const auto& pos = mesh->vertex(v)->pos;
        EXPECT_TRUE(std::isfinite(pos[0])) << "vertex " << v << " u not finite";
        EXPECT_TRUE(std::isfinite(pos[1])) << "vertex " << v << " v not finite";
        EXPECT_FLOAT_EQ(pos[2], 0.f) << "vertex " << v << " z != 0";
    }

    // Check no triangle flips (all faces have positive signed area in UV space)
    for (const auto& f : mesh->faces()) {
        auto e = f->head;
        const auto& p0 = e->vertex->pos;
        const auto& p1 = e->next->vertex->pos;
        const auto& p2 = e->next->next->vertex->pos;
        auto area = (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]);
        EXPECT_GT(area, 0.f) << "face " << f->idx << " is flipped";
    }
}

TEST(HLSCM, WavySurface)
{
    // Spatially varying curvature: exercises multi-level hierarchy
    using HLSCM = HierarchicalLSCM<float>;
    auto mesh = ConstructWavySurface<HLSCM::Mesh>(20, 20);
    HLSCM::Compute(mesh);

    for (std::size_t v = 0; v < mesh->num_vertices(); ++v) {
        const auto& pos = mesh->vertex(v)->pos;
        EXPECT_TRUE(std::isfinite(pos[0])) << "vertex " << v << " u not finite";
        EXPECT_TRUE(std::isfinite(pos[1])) << "vertex " << v << " v not finite";
        EXPECT_FLOAT_EQ(pos[2], 0.f) << "vertex " << v << " z != 0";
    }

    // Check no triangle flips
    for (const auto& f : mesh->faces()) {
        auto e = f->head;
        const auto& p0 = e->vertex->pos;
        const auto& p1 = e->next->vertex->pos;
        const auto& p2 = e->next->next->vertex->pos;
        auto area = (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]);
        EXPECT_GT(area, 0.f) << "face " << f->idx << " is flipped";
    }
}

TEST(HLSCM, InstanceAPILevelRatio)
{
    // Verify setLevelRatio/setMinCoarseVertices produce valid results
    using HLSCM = HierarchicalLSCM<float>;

    auto mesh = ConstructWavySurface<HLSCM::Mesh>(20, 20);
    HLSCM hlscm;
    hlscm.setLevelRatio(4);
    hlscm.setMinCoarseVertices(50);
    hlscm.compute(mesh);

    for (std::size_t v = 0; v < mesh->num_vertices(); ++v) {
        const auto& pos = mesh->vertex(v)->pos;
        EXPECT_TRUE(std::isfinite(pos[0])) << "vertex " << v;
        EXPECT_TRUE(std::isfinite(pos[1])) << "vertex " << v;
        EXPECT_FLOAT_EQ(pos[2], 0.f);
    }
}

TEST(HLSCM, MultiLevelHierarchy)
{
    // Verify that a sufficiently large mesh actually triggers multi-level hierarchy
    // by checking that HLSCM produces different results from single-level LSCM
    // when using ABF-optimized angles (the finest-level angle preservation matters)
    using ABFType = ABFPlusPlus<float>;
    using HLSCM = HierarchicalLSCM<float, ABFType::Mesh>;
    using LSCM = AngleBasedLSCM<float, ABFType::Mesh>;

    auto mesh_lscm = ConstructWavySurface<ABFType::Mesh>(20, 20);
    ABFType::Compute(mesh_lscm);
    LSCM::Compute(mesh_lscm);

    auto mesh_hlscm = ConstructWavySurface<ABFType::Mesh>(20, 20);
    ABFType::Compute(mesh_hlscm);
    HLSCM::Compute(mesh_hlscm);

    // Both should produce valid UVs
    for (std::size_t v = 0; v < mesh_hlscm->num_vertices(); ++v) {
        const auto& pos = mesh_hlscm->vertex(v)->pos;
        EXPECT_TRUE(std::isfinite(pos[0])) << "vertex " << v;
        EXPECT_TRUE(std::isfinite(pos[1])) << "vertex " << v;
        EXPECT_FLOAT_EQ(pos[2], 0.f);
    }

    // With a multi-level hierarchy, the cascadic solve path differs from
    // single-level LSCM, so results will differ slightly (both are valid)
    bool anyDifference = false;
    for (std::size_t v = 0; v < mesh_lscm->num_vertices(); ++v) {
        for (int i = 0; i < 2; ++i) {
            if (std::abs(mesh_hlscm->vertex(v)->pos[i] - mesh_lscm->vertex(v)->pos[i]) > 1e-4f) {
                anyDifference = true;
                break;
            }
        }
        if (anyDifference)
            break;
    }
    EXPECT_TRUE(anyDifference) << "HLSCM and LSCM produced identical results — "
                                  "multi-level hierarchy may not have been triggered";
}

TEST(HLSCM, ABFPlusPlusAnglePreservation)
{
    // Verify that ABF++ → HLSCM uses the ABF-optimized angles at the finest
    // hierarchy level by checking that the result differs from geometry-only
    // HLSCM.  If angles were being discarded (recomputed from geometry), both
    // would produce identical output.
    using ABFType = ABFPlusPlus<float>;
    using HLSCM_ABF = HierarchicalLSCM<float, ABFType::Mesh>;
    using HLSCM_Geo = HierarchicalLSCM<float>;

    constexpr std::size_t N = 20;

    // HLSCM with ABF-optimized angles
    auto mesh_abf = ConstructWavySurface<ABFType::Mesh>(N, N);
    ABFType::Compute(mesh_abf);
    HLSCM_ABF::Compute(mesh_abf);

    // HLSCM with geometry-only angles (no ABF)
    auto mesh_geo = ConstructWavySurface<HLSCM_Geo::Mesh>(N, N);
    HLSCM_Geo::Compute(mesh_geo);

    // Both must produce valid UVs
    for (std::size_t v = 0; v < mesh_abf->num_vertices(); ++v) {
        const auto& pos = mesh_abf->vertex(v)->pos;
        EXPECT_TRUE(std::isfinite(pos[0])) << "vertex " << v;
        EXPECT_TRUE(std::isfinite(pos[1])) << "vertex " << v;
        EXPECT_FLOAT_EQ(pos[2], 0.f);
    }

    // The ABF-optimized angles should produce a measurably different result
    // from geometry angles, proving they are actually being used at the
    // finest level.
    bool anyDifference = false;
    for (std::size_t v = 0; v < mesh_abf->num_vertices(); ++v) {
        for (int i = 0; i < 2; ++i) {
            if (std::abs(mesh_abf->vertex(v)->pos[i] - mesh_geo->vertex(v)->pos[i]) > 1e-6f) {
                anyDifference = true;
                break;
            }
        }
        if (anyDifference)
            break;
    }
    EXPECT_TRUE(anyDifference) << "ABF++ angles had no effect on HLSCM output — "
                                  "angles may not be preserved at finest hierarchy level";
}

TEST(HLSCM, ABFReducesConformalDistortion)
{
    // ABF-optimized angles should always reduce conformal distortion in LSCM,
    // including HLSCM.  Measure angle error between the input 3D mesh and the
    // flattened 2D parameterization: for each edge, compare the original 3D
    // interior angle with the corresponding 2D angle.  ABF++ → HLSCM should
    // have lower total angle error than geometry-only HLSCM.
    using ABFType = ABFPlusPlus<float>;
    using HLSCM_ABF = HierarchicalLSCM<float, ABFType::Mesh>;
    using HLSCM_Geo = HierarchicalLSCM<float>;

    constexpr std::size_t N = 20;

    // Helper: save per-edge 3D angles keyed by edge idx.
    // Note: num_edges() counts only face-adjacent half-edges, but edge idx
    // values are assigned from the full half-edge pool (including boundary
    // half-edges), so they are not contiguous.  Use a map to avoid OOB.
    auto save3DAngles = [](const auto& mesh) {
        std::unordered_map<std::size_t, float> angles;
        for (const auto& f : mesh->faces()) {
            for (auto& e : *f) {
                angles[e->idx] = e->alpha;
            }
        }
        return angles;
    };

    // Helper: compute total angle distortion (sum of squared angle errors)
    auto angleDistortion = [](const auto& mesh,
                              const std::unordered_map<std::size_t, float>& origAngles) {
        double totalErr = 0.0;
        for (const auto& f : mesh->faces()) {
            for (auto& e : *f) {
                auto ab = e->next->vertex->pos - e->vertex->pos;
                auto ac = e->next->next->vertex->pos - e->vertex->pos;
                auto uvAngle = OpenABF::interior_angle(ab, ac);
                double diff = static_cast<double>(uvAngle) - origAngles.at(e->idx);
                totalErr += diff * diff;
            }
        }
        return totalErr;
    };

    // ABF++ → HLSCM path
    auto mesh_abf = ConstructWavySurface<ABFType::Mesh>(N, N);
    ABFType::Compute(mesh_abf);
    auto angles_abf = save3DAngles(mesh_abf);
    HLSCM_ABF::Compute(mesh_abf);
    auto distortion_abf = angleDistortion(mesh_abf, angles_abf);

    // Geometry-only HLSCM path
    auto mesh_geo = ConstructWavySurface<HLSCM_Geo::Mesh>(N, N);
    auto angles_geo = save3DAngles(mesh_geo);
    HLSCM_Geo::Compute(mesh_geo);
    auto distortion_geo = angleDistortion(mesh_geo, angles_geo);

    EXPECT_LT(distortion_abf, distortion_geo)
        << "ABF++ did not reduce conformal distortion in HLSCM: "
        << "ABF=" << distortion_abf << " vs Geo=" << distortion_geo;
}

TEST(HLSCM, PerformanceComparison)
{
    // Compare HLSCM vs AngleBasedLSCM (with same LSCG solver) on a large mesh.
    // Both use LeastSquaresConjugateGradient so the comparison isolates the
    // benefit of the hierarchical initial guess.
    using SolverType = Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<float>>;
    using HLSCM = HierarchicalLSCM<float>;
    using LSCM = AngleBasedLSCM<float, HalfEdgeMesh<float>, SolverType>;

    constexpr std::size_t rows = 75;
    constexpr std::size_t cols = 75;
    // 5625 vertices, 10952 faces — large enough for meaningful hierarchy

    // Time AngleBasedLSCM (with LSCG solver, no hierarchy)
    auto mesh_lscm = ConstructWavySurface<LSCM::Mesh>(rows, cols);
    auto t0 = std::chrono::steady_clock::now();
    LSCM::Compute(mesh_lscm);
    auto t1 = std::chrono::steady_clock::now();
    auto lscm_us = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();

    // Time HierarchicalLSCM
    auto mesh_hlscm = ConstructWavySurface<HLSCM::Mesh>(rows, cols);
    auto t2 = std::chrono::steady_clock::now();
    HLSCM::Compute(mesh_hlscm);
    auto t3 = std::chrono::steady_clock::now();
    auto hlscm_us = std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count();

    std::cout << "\n  HLSCM Performance (" << rows << "x" << cols << " wavy surface, "
              << mesh_hlscm->num_vertices() << " verts, " << mesh_hlscm->num_faces() << " faces):\n"
              << "    LSCM (LSCG, no hierarchy): " << lscm_us / 1000.0 << " ms\n"
              << "    HLSCM (hierarchical):      " << hlscm_us / 1000.0 << " ms\n"
              << "    Speedup:                   "
              << static_cast<double>(lscm_us) /
                     static_cast<double>(std::max(hlscm_us, decltype(hlscm_us)(1)))
              << "x\n";

    // Verify HLSCM produced valid results
    for (std::size_t v = 0; v < mesh_hlscm->num_vertices(); ++v) {
        const auto& pos = mesh_hlscm->vertex(v)->pos;
        EXPECT_TRUE(std::isfinite(pos[0])) << "vertex " << v;
        EXPECT_TRUE(std::isfinite(pos[1])) << "vertex " << v;
        EXPECT_FLOAT_EQ(pos[2], 0.f) << "vertex " << v;
    }

    // Verify no triangle flips in HLSCM output
    for (const auto& f : mesh_hlscm->faces()) {
        auto e = f->head;
        const auto& p0 = e->vertex->pos;
        const auto& p1 = e->next->vertex->pos;
        const auto& p2 = e->next->next->vertex->pos;
        auto area = (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]);
        EXPECT_GT(area, 0.f) << "face " << f->idx << " is flipped";
    }
}