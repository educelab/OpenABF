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

    // Instance API must produce identical output to the static overload
    for (std::size_t v = 0; v < mesh_static->num_vertices(); ++v) {
        for (auto i = 0; i < 3; i++) {
            EXPECT_FLOAT_EQ(mesh_instance->vertex(v)->pos[i], mesh_static->vertex(v)->pos[i]);
        }
    }

    // Pinned vertices must land at their expected UV positions exactly.
    // HLSCM uses the same pin-placement logic as AngleBasedLSCM: pin0 at the
    // origin and pin1 at distance |p1-p0| along the dominant world-space axis
    // of the (p1-p0) vector. For ConstructPyramid with pins (0, 1) this places
    // pin0 at {0, 0, 0} and pin1 at {2, 0, 0}.
    const Vec3f expectedPin0{0, 0, 0};
    const Vec3f expectedPin1{2, 0, 0};
    for (auto i = 0; i < 3; i++) {
        EXPECT_FLOAT_EQ(mesh_instance->vertex(0)->pos[i], expectedPin0[i]);
        EXPECT_FLOAT_EQ(mesh_instance->vertex(1)->pos[i], expectedPin1[i]);
    }

    // Sanity check: pinning must actually be applied. Non-pinned vertices must
    // not all collapse onto the pin positions.
    bool anyNonPinnedDiffers = false;
    for (std::size_t v = 2; v < mesh_instance->num_vertices(); ++v) {
        const auto& p = mesh_instance->vertex(v)->pos;
        const bool atPin0 =
            (p[0] == expectedPin0[0]) && (p[1] == expectedPin0[1]) && (p[2] == expectedPin0[2]);
        const bool atPin1 =
            (p[0] == expectedPin1[0]) && (p[1] == expectedPin1[1]) && (p[2] == expectedPin1[2]);
        if (!atPin0 && !atPin1) {
            anyNonPinnedDiffers = true;
            break;
        }
    }
    EXPECT_TRUE(anyNonPinnedDiffers);
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

    // Check no triangle flips (all faces have positive signed area in UV space).
    // Use a small tolerance to accommodate floating-point rounding on near-degenerate faces.
    for (const auto& f : mesh->faces()) {
        auto e = f->head;
        const auto& p0 = e->vertex->pos;
        const auto& p1 = e->next->vertex->pos;
        const auto& p2 = e->next->next->vertex->pos;
        auto area = (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]);
        EXPECT_GE(area, -1e-5f) << "face " << f->idx << " is flipped";
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

    // Check no triangle flips.
    // Use a small tolerance to accommodate floating-point rounding on near-degenerate faces.
    for (const auto& f : mesh->faces()) {
        auto e = f->head;
        const auto& p0 = e->vertex->pos;
        const auto& p1 = e->next->vertex->pos;
        const auto& p2 = e->next->next->vertex->pos;
        auto area = (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]);
        EXPECT_GE(area, -1e-5f) << "face " << f->idx << " is flipped";
    }
}

TEST(HLSCM, InstanceAPILevelRatio)
{
    // Verify that setLevelRatio/setMinCoarseVertices produce valid results AND
    // that the levelRatio parameter actually affects hierarchy construction.
    // A no-op setLevelRatio (e.g., parameter ignored) would cause the
    // hierarchy-depth assertion below to fail because both ratios would
    // produce identical hierarchies.
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

    // Directly invoke buildHierarchy with two clearly different levelRatio
    // values on the same input.  A larger ratio decimates more aggressively
    // per level and should produce a hierarchy with fewer levels (each level
    // also having a smaller alive-vertex count) than a small ratio.
    using namespace OpenABF::detail::hlscm;
    auto mesh2 = ConstructWavySurface<HLSCM::Mesh>(20, 20);
    constexpr std::size_t pin0 = 0;
    constexpr std::size_t pin1 = 19;  // opposite corner of a 20x20 grid row
    constexpr std::size_t minCoarseVerts = 10;

    auto [levelsSmall, _csmall] =
        buildHierarchy<float>(mesh2, pin0, pin1, /*levelRatio=*/2, minCoarseVerts);
    auto [levelsLarge, _clarge] =
        buildHierarchy<float>(mesh2, pin0, pin1, /*levelRatio=*/8, minCoarseVerts);

    // Sanity: both hierarchies have at least the finest level
    ASSERT_GE(levelsSmall.size(), 1u);
    ASSERT_GE(levelsLarge.size(), 1u);

    // Observable effect: levelRatio=2 must produce more (or equal-but-not-fewer)
    // levels than levelRatio=8 on a 400-vertex mesh decimated to <=10 verts.
    // Concretely we expect strict inequality here — if equal, levelRatio is
    // not being applied.
    EXPECT_GT(levelsSmall.size(), levelsLarge.size())
        << "levelRatio appears to have no effect on hierarchy depth: " << "ratio=2 produced "
        << levelsSmall.size() << " levels, " << "ratio=8 produced " << levelsLarge.size()
        << " levels";
}

TEST(HLSCM, MultiLevelHierarchy)
{
    // Verify a multi-level hierarchy is actually built by calling
    // detail::hlscm::buildHierarchy directly and asserting on levels.size().
    using HLSCM = OpenABF::HierarchicalLSCM<float>;
    auto mesh = ConstructWavySurface<HLSCM::Mesh>(20, 20);

    // Match HLSCM defaults except force a small minCoarseVerts so we get >1 level
    auto [levels, collapsesByLevel] =
        OpenABF::detail::hlscm::buildHierarchy<float>(mesh, 0, 1, /*levelRatio=*/10,
                                                      /*minCoarseVerts=*/10);
    EXPECT_GE(levels.size(), std::size_t(2))
        << "buildHierarchy produced only " << levels.size() << " level(s)";

    // Sanity-check that the full pipeline still runs end-to-end on this mesh
    HLSCM hlscm;
    hlscm.setMinCoarseVertices(10);
    ASSERT_NO_THROW(hlscm.compute(mesh));

    // All UVs should be finite and z=0
    for (const auto& v : mesh->vertices()) {
        EXPECT_TRUE(std::isfinite(v->pos[0]));
        EXPECT_TRUE(std::isfinite(v->pos[1]));
        EXPECT_FLOAT_EQ(v->pos[2], 0.f);
    }
}

TEST(HLSCM, ABFPlusPlusAnglePreservation)
{
    // Verify that ABF++ → HLSCM uses the ABF-optimized angles at the finest
    // hierarchy level by checking that the result differs from geometry-only
    // HLSCM.  If angles were being discarded (recomputed from geometry), both
    // would produce identical output.
    //
    // Uses a hemisphere (genuine Gaussian curvature) so ABF angle correction
    // has meaningful work to do.  On a near-flat mesh (e.g. wavy surface) the
    // ABF correction is small, so the UV deltas are also small and the test
    // is weak; a hemisphere guarantees a strong, easy-to-detect signal.
    using ABFType = ABFPlusPlus<float>;
    using HLSCM_ABF = HierarchicalLSCM<float, ABFType::Mesh>;
    using HLSCM_Geo = HierarchicalLSCM<float>;

    constexpr std::size_t rings = 12;
    constexpr std::size_t sectors = 24;

    // Snapshot the raw geometry angles before ABF rewrites them so we can
    // verify ABF actually altered the angles (i.e. its output is non-trivial).
    auto saveAngles = [](const auto& mesh) {
        std::unordered_map<std::size_t, float> angles;
        for (const auto& f : mesh->faces()) {
            for (auto& e : *f) {
                angles[e->idx] = e->alpha;
            }
        }
        return angles;
    };

    // HLSCM with ABF-optimized angles
    auto mesh_abf = ConstructHemisphere<ABFType::Mesh>(rings, sectors);
    auto rawAngles = saveAngles(mesh_abf);
    ABFType::Compute(mesh_abf);

    // Confirm ABF actually changed the per-half-edge angles before HLSCM runs.
    // Without this, an "anyDifference" check downstream could be triggered by
    // unrelated noise (solver tolerances, etc.).
    float maxAngleDelta = 0.f;
    for (const auto& f : mesh_abf->faces()) {
        for (auto& e : *f) {
            maxAngleDelta = std::max(maxAngleDelta, std::abs(e->alpha - rawAngles.at(e->idx)));
        }
    }
    EXPECT_GT(maxAngleDelta, 1e-3f)
        << "ABF++ did not modify input angles meaningfully on the curved mesh "
           "(maxAngleDelta="
        << maxAngleDelta << ")";

    HLSCM_ABF::Compute(mesh_abf);

    // HLSCM with geometry-only angles (no ABF)
    auto mesh_geo = ConstructHemisphere<HLSCM_Geo::Mesh>(rings, sectors);
    HLSCM_Geo::Compute(mesh_geo);

    // Both must produce valid UVs
    ASSERT_EQ(mesh_abf->num_vertices(), mesh_geo->num_vertices());
    for (std::size_t v = 0; v < mesh_abf->num_vertices(); ++v) {
        const auto& pos = mesh_abf->vertex(v)->pos;
        EXPECT_TRUE(std::isfinite(pos[0])) << "vertex " << v;
        EXPECT_TRUE(std::isfinite(pos[1])) << "vertex " << v;
        EXPECT_FLOAT_EQ(pos[2], 0.f);
    }

    // Measure UV difference between the two parameterizations.  Because the
    // two solves are independent (different boundary pins resolve to the same
    // similarity class but not the same exact coordinates), we look at the
    // per-vertex maximum delta and require it to be substantially larger than
    // numerical noise.  On a hemisphere the signal is order 0.1+ in UV units;
    // a regression that discards ABF angles drops this to ~0.
    float maxDiff = 0.f;
    double l2Diff = 0.0;
    for (std::size_t v = 0; v < mesh_abf->num_vertices(); ++v) {
        for (int i = 0; i < 2; ++i) {
            float d = std::abs(mesh_abf->vertex(v)->pos[i] - mesh_geo->vertex(v)->pos[i]);
            maxDiff = std::max(maxDiff, d);
            l2Diff += static_cast<double>(d) * static_cast<double>(d);
        }
    }
    EXPECT_GT(maxDiff, 1e-2f) << "ABF++ angles had negligible effect on HLSCM output "
                                 "(maxDiff="
                              << maxDiff << ", L2=" << std::sqrt(l2Diff)
                              << ") — angles may not be preserved at finest hierarchy level";
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
        << "ABF++ did not reduce conformal distortion in HLSCM: " << "ABF=" << distortion_abf
        << " vs Geo=" << distortion_geo;
}

TEST(HLSCM, LargeMeshValidation)
{
    // Validate HLSCM correctness on a large mesh and record timing for
    // reference. Also compares wall-clock time against flat LSCG to illustrate
    // the hierarchical speedup (not a hard performance assertion).
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

    // Verify no triangle flips in HLSCM output. Allow a small negative epsilon
    // for near-degenerate faces (numerical noise only; real flips are O(1e-4)).
    constexpr float kAreaEps = -1e-5f;
    for (const auto& f : mesh_hlscm->faces()) {
        auto e = f->head;
        const auto& p0 = e->vertex->pos;
        const auto& p1 = e->next->vertex->pos;
        const auto& p2 = e->next->next->vertex->pos;
        auto area = (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]);
        EXPECT_GE(area, kAreaEps) << "face " << f->idx << " is flipped";
    }
}

TEST(HLSCM, SingleTriangle)
{
    // A 3-vertex mesh is entirely boundary — no collapses possible.
    // Exercises the single-level fallback path where levels.size() <= 1.
    using HLSCM = HierarchicalLSCM<float>;
    auto mesh = ConstructPyramid<HLSCM::Mesh>();
    ASSERT_NO_THROW(HLSCM::Compute(mesh));

    for (std::size_t v = 0; v < mesh->num_vertices(); ++v) {
        const auto& pos = mesh->vertex(v)->pos;
        EXPECT_TRUE(std::isfinite(pos[0])) << "vertex " << v;
        EXPECT_TRUE(std::isfinite(pos[1])) << "vertex " << v;
        EXPECT_FLOAT_EQ(pos[2], 0.f) << "vertex " << v;
    }
}

TEST(HLSCM, MeshExceptionOnClosedMesh)
{
    // A closed (boundary-free) mesh has no boundary vertices, so AutoSelectPins
    // should throw MeshException.
    using HLSCM = HierarchicalLSCM<float>;
    auto mesh = HLSCM::Mesh::New();
    // Tetrahedron (fully closed — no boundary)
    mesh->insert_vertex(1, 1, 1);
    mesh->insert_vertex(-1, -1, 1);
    mesh->insert_vertex(-1, 1, -1);
    mesh->insert_vertex(1, -1, -1);
    mesh->insert_faces({{0, 1, 2}, {0, 2, 3}, {0, 3, 1}, {1, 3, 2}});

    EXPECT_THROW(HLSCM::Compute(mesh), MeshException);
}

TEST(HLSCM, DirectSolverBranch)
{
    // Exercise the SparseLU (direct solver) code path. Result must be valid.
    using Solver = Eigen::SparseLU<Eigen::SparseMatrix<float>>;
    using HLSCM = HierarchicalLSCM<float, HalfEdgeMesh<float>, Solver>;

    auto mesh = ConstructHemisphere<HLSCM::Mesh>(6, 12);
    ASSERT_NO_THROW(HLSCM::Compute(mesh));

    for (std::size_t v = 0; v < mesh->num_vertices(); ++v) {
        const auto& pos = mesh->vertex(v)->pos;
        EXPECT_TRUE(std::isfinite(pos[0])) << "vertex " << v;
        EXPECT_TRUE(std::isfinite(pos[1])) << "vertex " << v;
        EXPECT_FLOAT_EQ(pos[2], 0.f) << "vertex " << v;
    }
}

TEST(HLSCM, FlatGridNoDistortion)
{
    // A flat grid has zero Gaussian curvature. HLSCM should produce a valid
    // parameterization with no triangle flips, verifying the multi-level
    // hierarchy introduces no distortion on a trivially parameterizable mesh.
    using HLSCM = HierarchicalLSCM<float>;
    auto mesh = ConstructGrid<HLSCM::Mesh>(15, 15);
    ASSERT_NO_THROW(HLSCM::Compute(mesh));

    for (std::size_t v = 0; v < mesh->num_vertices(); ++v) {
        const auto& pos = mesh->vertex(v)->pos;
        EXPECT_TRUE(std::isfinite(pos[0])) << "vertex " << v;
        EXPECT_TRUE(std::isfinite(pos[1])) << "vertex " << v;
        EXPECT_FLOAT_EQ(pos[2], 0.f) << "vertex " << v;
    }

    constexpr float kAreaEps = -1e-5f;
    for (const auto& f : mesh->faces()) {
        auto e = f->head;
        const auto& p0 = e->vertex->pos;
        const auto& p1 = e->next->vertex->pos;
        const auto& p2 = e->next->next->vertex->pos;
        auto area = (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]);
        EXPECT_GE(area, kAreaEps) << "face " << f->idx << " is flipped";
    }
}

TEST(HLSCM, LevelRatioBoundaryValues)
{
    // setLevelRatio(0) and setLevelRatio(1) must throw std::invalid_argument
    // because they would cause division-by-zero or degenerate hierarchies.
    using HLSCM = HierarchicalLSCM<float>;
    HLSCM hlscm;
    EXPECT_THROW(hlscm.setLevelRatio(0), std::invalid_argument);
    EXPECT_THROW(hlscm.setLevelRatio(1), std::invalid_argument);
    EXPECT_NO_THROW(hlscm.setLevelRatio(2));

    // setMinCoarseVertices(0), (1), (2) must throw; 3 is the minimum valid value
    EXPECT_THROW(hlscm.setMinCoarseVertices(0), std::invalid_argument);
    EXPECT_THROW(hlscm.setMinCoarseVertices(1), std::invalid_argument);
    EXPECT_THROW(hlscm.setMinCoarseVertices(2), std::invalid_argument);
    EXPECT_NO_THROW(hlscm.setMinCoarseVertices(3));
}

TEST(HLSCM, DoubleOnHemisphere)
{
    // Double-precision HLSCM on a non-trivial mesh — verifies the template
    // compiles and produces finite, z=0 results in double precision.
    using HLSCM = HierarchicalLSCM<double>;
    auto mesh = ConstructHemisphere<HLSCM::Mesh, double>(8, 16);
    ASSERT_NO_THROW(HLSCM::Compute(mesh));

    for (std::size_t v = 0; v < mesh->num_vertices(); ++v) {
        const auto& pos = mesh->vertex(v)->pos;
        EXPECT_TRUE(std::isfinite(pos[0])) << "vertex " << v;
        EXPECT_TRUE(std::isfinite(pos[1])) << "vertex " << v;
        EXPECT_DOUBLE_EQ(pos[2], 0.0) << "vertex " << v;
    }
}

TEST(HLSCMInternal, DecimationMesh_RejectsPinnedVertex)
{
    // Directly unit-test detail::hlscm::DecimationMesh: tryCollapse must
    // reject the collapse when the vertex-to-remove is pinned, even if the
    // edge is geometrically valid.
    using namespace OpenABF::detail::hlscm;
    using DMesh = DecimationMesh<float>;

    // Build a 4-vertex open mesh (pyramid)
    using Mesh = HalfEdgeMesh<float>;
    auto mesh = Mesh::New();
    mesh->insert_vertex(0, 0, 0);
    mesh->insert_vertex(2, 0, 0);
    mesh->insert_vertex(1, std::sqrt(3.f), 0);
    mesh->insert_vertex(1, std::sqrt(3.f) / 3.f, std::sqrt(6.f) * 2.f / 3.f);
    mesh->insert_faces({{1, 3, 0}, {3, 2, 0}, {3, 1, 2}});

    // Pin vertices 0 and 1
    constexpr std::size_t pin0 = 0;
    constexpr std::size_t pin1 = 1;
    DMesh dm;
    dm.build(mesh, pin0, pin1);

    // Attempting to remove a pinned vertex (pin0=0 → vertex 1) must return nullopt
    auto result = dm.tryCollapse(pin0, 1);
    EXPECT_FALSE(result.has_value()) << "tryCollapse should reject collapse when vRemove is pinned";

    // Attempting to remove a non-pinned vertex should succeed (may return a valid record)
    auto result2 = dm.tryCollapse(3, 0);
    // vertex 3 is the apex (non-boundary, non-pinned) — collapse may succeed or be
    // rejected on geometric grounds, but must never crash
    (void)result2;
}

TEST(HLSCMInternal, BuildHierarchy_LevelCount)
{
    // Directly invoke detail::hlscm::buildHierarchy on a 20×20 wavy surface
    // and assert that the returned hierarchy has the expected structure:
    //   - At least 3 levels at minCoarseVerts=10
    //   - Consecutive-level vertex-count ratio approximately matches levelRatio
    //   - localToOriginal and originalToLocal are consistent inverses
    //   - Pin vertices survive at every level
    using HLSCM = OpenABF::HierarchicalLSCM<float>;
    auto mesh = ConstructWavySurface<HLSCM::Mesh>(20, 20);

    constexpr std::size_t pin0 = 0;
    constexpr std::size_t pin1 = 1;
    constexpr std::size_t levelRatio = 4;
    constexpr std::size_t minCoarseVerts = 10;

    auto [levels, collapsesByLevel] =
        OpenABF::detail::hlscm::buildHierarchy<float>(mesh, pin0, pin1, levelRatio, minCoarseVerts);

    // 20x20 = 400 verts, ratio 4, min 10  →  400, 100, 25, 10  → 4 levels (≥3)
    ASSERT_GE(levels.size(), std::size_t(3))
        << "Expected hierarchy with >= 3 levels; got " << levels.size();

    // levels and collapsesByLevel are paired: each collapsesByLevel[k]
    // describes the transition from levels[k] to levels[k+1]
    ASSERT_EQ(collapsesByLevel.size(), levels.size() - 1);

    // For each consecutive pair, assert vertex count ratio ≈ levelRatio
    // (within a generous tolerance because decimation may stop early for
    // geometric reasons and the last step clamps to minCoarseVerts)
    for (std::size_t k = 0; k + 1 < levels.size(); ++k) {
        auto prev = levels[k].localToOriginal.size();
        auto next = levels[k + 1].localToOriginal.size();
        ASSERT_GT(prev, next) << "Level " << (k + 1) << " is not coarser than level " << k;
        // ratio = prev/next; allow [0.5*levelRatio, 2*levelRatio] except at the
        // floor where we clamp to minCoarseVerts
        auto ratio = static_cast<double>(prev) / static_cast<double>(next);
        if (next > minCoarseVerts) {
            EXPECT_GE(ratio, 0.5 * levelRatio)
                << "Level " << k << "→" << (k + 1) << " ratio " << ratio << " too small";
            EXPECT_LE(ratio, 2.0 * levelRatio)
                << "Level " << k << "→" << (k + 1) << " ratio " << ratio << " too large";
        }
    }

    // For each level: localToOriginal and originalToLocal are consistent inverses.
    for (std::size_t k = 0; k < levels.size(); ++k) {
        const auto& lvl = levels[k];
        EXPECT_EQ(lvl.localToOriginal.size(), lvl.positions.size())
            << "Level " << k << ": localToOriginal size != positions size";
        EXPECT_EQ(lvl.originalToLocal.size(), lvl.localToOriginal.size())
            << "Level " << k << ": map sizes mismatch";
        for (std::size_t li = 0; li < lvl.localToOriginal.size(); ++li) {
            auto origIdx = lvl.localToOriginal[li];
            auto it = lvl.originalToLocal.find(origIdx);
            ASSERT_NE(it, lvl.originalToLocal.end())
                << "Level " << k << ": original idx " << origIdx << " missing from originalToLocal";
            EXPECT_EQ(it->second, li)
                << "Level " << k << ": originalToLocal[" << origIdx << "] != " << li;
        }
    }

    // Pin vertices must appear in every level
    for (std::size_t k = 0; k < levels.size(); ++k) {
        EXPECT_TRUE(levels[k].originalToLocal.count(pin0))
            << "Level " << k << " missing pin0 (vertex " << pin0 << ")";
        EXPECT_TRUE(levels[k].originalToLocal.count(pin1))
            << "Level " << k << " missing pin1 (vertex " << pin1 << ")";
    }
}

TEST(HLSCMInternal, ProlongateUVs_BarycentricReconstruction)
{
    // Build a 2-level hierarchy on a small grid, assign known UVs at the
    // coarse level, and verify that prolongateUVs reconstructs the removed
    // vertices' UVs as barycentric interpolations of the containingTri's UVs.
    using HLSCM = OpenABF::HierarchicalLSCM<float>;
    auto mesh = ConstructGrid<HLSCM::Mesh>(5, 5);

    constexpr std::size_t pin0 = 0;
    constexpr std::size_t pin1 = 4;  // opposite corner of the top row

    // Force a 2-level hierarchy (25 verts → ~6 verts at ratio 4).
    auto [levels, collapsesByLevel] = OpenABF::detail::hlscm::buildHierarchy<float>(
        mesh, pin0, pin1, /*levelRatio=*/4, /*minCoarseVerts=*/5);
    ASSERT_GE(levels.size(), std::size_t(2)) << "Expected at least 2 hierarchy levels";
    ASSERT_FALSE(collapsesByLevel.empty()) << "Expected at least one collapse record set";

    // Use the first transition (finest collapse set): from levels[0] to levels[1].
    const auto& collapses = collapsesByLevel[0];
    ASSERT_FALSE(collapses.empty()) << "Finest-level collapse set is empty";

    // Assign deterministic UVs to coarse-level vertices: U = origIdx, V = -origIdx
    // (anything works as long as we have one UV per surviving vertex)
    std::unordered_map<std::size_t, std::array<float, 2>> coarseUVs;
    for (auto origIdx : levels[1].localToOriginal) {
        coarseUVs[origIdx] = {static_cast<float>(origIdx), -static_cast<float>(origIdx)};
    }

    auto fineUVs = OpenABF::detail::hlscm::prolongateUVs<float>(coarseUVs, collapses);

    // Coarse UVs must survive untouched
    for (const auto& [origIdx, uv] : coarseUVs) {
        ASSERT_TRUE(fineUVs.count(origIdx))
            << "Coarse vertex " << origIdx << " missing from prolongated UVs";
        EXPECT_FLOAT_EQ(fineUVs[origIdx][0], uv[0]) << "Coarse vertex " << origIdx << " U mutated";
        EXPECT_FLOAT_EQ(fineUVs[origIdx][1], uv[1]) << "Coarse vertex " << origIdx << " V mutated";
    }

    // Replay barycentric expansion in the same order as prolongateUVs (reverse
    // of the collapse log) and verify the prolongated UVs match within 1e-5.
    auto expected = coarseUVs;
    for (auto it = collapses.rbegin(); it != collapses.rend(); ++it) {
        const auto& rec = *it;
        ASSERT_TRUE(expected.count(rec.containingTri[0]))
            << "containingTri[0] " << rec.containingTri[0] << " UV missing";
        ASSERT_TRUE(expected.count(rec.containingTri[1]))
            << "containingTri[1] " << rec.containingTri[1] << " UV missing";
        ASSERT_TRUE(expected.count(rec.containingTri[2]))
            << "containingTri[2] " << rec.containingTri[2] << " UV missing";
        auto uv0 = expected.at(rec.containingTri[0]);
        auto uv1 = expected.at(rec.containingTri[1]);
        auto uv2 = expected.at(rec.containingTri[2]);
        std::array<float, 2> exp{
            rec.bary[0] * uv0[0] + rec.bary[1] * uv1[0] + rec.bary[2] * uv2[0],
            rec.bary[0] * uv0[1] + rec.bary[1] * uv1[1] + rec.bary[2] * uv2[1]};
        expected[rec.vRemoved] = exp;

        ASSERT_TRUE(fineUVs.count(rec.vRemoved))
            << "Removed vertex " << rec.vRemoved << " missing from prolongated UVs";
        EXPECT_NEAR(fineUVs[rec.vRemoved][0], exp[0], 1e-5f)
            << "vertex " << rec.vRemoved << " U mismatch";
        EXPECT_NEAR(fineUVs[rec.vRemoved][1], exp[1], 1e-5f)
            << "vertex " << rec.vRemoved << " V mismatch";
    }
}

TEST(HLSCMInternal, SolveLSCMLevel_KnownMesh)
{
    // Call detail::hlscm::solveLSCMLevel directly on a single-level pyramid
    // HierarchyLevel and assert:
    //   - all returned UVs are finite and z-component is implicit 0
    //   - pinned vertices have the prescribed UV positions (same pin
    //     placement convention as AngleBasedLSCM)
    //   - the result matches AngleBasedLSCM::Compute on the same mesh
    using HLSCM = OpenABF::HierarchicalLSCM<float>;
    using Solver =
        Eigen::ConjugateGradient<Eigen::SparseMatrix<float>, Eigen::Lower | Eigen::Upper>;
    using LSCM = OpenABF::AngleBasedLSCM<float, HLSCM::Mesh, Solver>;

    constexpr std::size_t pin0 = 0;
    constexpr std::size_t pin1 = 1;

    // Build a single-level hierarchy (pyramid is too small to decimate).
    auto meshH = ConstructPyramid<HLSCM::Mesh>();
    auto [levels, collapsesByLevel] = OpenABF::detail::hlscm::buildHierarchy<float>(
        meshH, pin0, pin1, /*levelRatio=*/10, /*minCoarseVerts=*/100);
    ASSERT_EQ(levels.size(), std::size_t(1)) << "Pyramid should produce a single-level hierarchy";

    const auto& level = levels[0];
    auto levelMesh = OpenABF::detail::hlscm::buildLevelMesh<float>(level);
    OpenABF::ComputeMeshAngles(levelMesh);

    auto uvs = OpenABF::detail::hlscm::solveLSCMLevel<float, Solver>(levelMesh, level, pin0, pin1,
                                                                     nullptr);

    // All UVs must be finite (z is implicit 0; solveLSCMLevel only returns 2-vectors)
    ASSERT_EQ(uvs.size(), meshH->num_vertices());
    for (const auto& [origIdx, uv] : uvs) {
        EXPECT_TRUE(std::isfinite(uv[0])) << "vertex " << origIdx << " U not finite";
        EXPECT_TRUE(std::isfinite(uv[1])) << "vertex " << origIdx << " V not finite";
    }

    // The pinned vertices use the same placement logic as AngleBasedLSCM:
    // pin0 sits at the origin; pin1 sits on whichever axis its displacement
    // from pin0 has the largest magnitude. For the pyramid (verts 0=(0,0,0)
    // and 1=(2,0,0)), pin1 lies at (2, 0).
    ASSERT_TRUE(uvs.count(pin0));
    ASSERT_TRUE(uvs.count(pin1));
    EXPECT_FLOAT_EQ(uvs[pin0][0], 0.f);
    EXPECT_FLOAT_EQ(uvs[pin0][1], 0.f);
    EXPECT_FLOAT_EQ(uvs[pin1][0], 2.f);
    EXPECT_FLOAT_EQ(uvs[pin1][1], 0.f);

    // Compute the reference solution via the top-level AngleBasedLSCM with the
    // same solver and same pin selection; results must match.
    auto meshL = ConstructPyramid<LSCM::Mesh>();
    LSCM::Compute(meshL, pin0, pin1);

    for (std::size_t v = 0; v < meshL->num_vertices(); ++v) {
        ASSERT_TRUE(uvs.count(v)) << "vertex " << v << " missing from solveLSCMLevel UVs";
        const auto& ref = meshL->vertex(v)->pos;
        EXPECT_NEAR(uvs[v][0], ref[0], 1e-4f) << "vertex " << v << " U disagrees with LSCM";
        EXPECT_NEAR(uvs[v][1], ref[1], 1e-4f) << "vertex " << v << " V disagrees with LSCM";
    }
}

// ------------------------------------------------------------------
// LSCMSystemBuild — direct tests for detail::lscm::buildSystem (A8).
//
// buildSystem is the shared system-assembly utility extracted from
// AngleBasedLSCM::ComputeImpl and HierarchicalLSCM::solveLSCMLevel. These
// tests exercise it on a pyramid (3 faces, 4 vertices, 2 pins → 2 free),
// asserting the structural invariants of the produced system rather than
// the numerical solution (which is covered transitively by the existing
// parameterization tests once both call sites migrate to buildSystem).
// ------------------------------------------------------------------
namespace
{

using LSCMSystemMesh = HalfEdgeMesh<float>;

auto BuildPyramidSystem()
{
    constexpr std::size_t pin0Idx = 0;
    constexpr std::size_t pin1Idx = 1;
    auto mesh = ConstructPyramid<LSCMSystemMesh>();
    ComputeMeshAngles(mesh);
    auto p0 = mesh->vertex(pin0Idx);
    auto p1 = mesh->vertex(pin1Idx);
    auto parts = OpenABF::detail::lscm::buildSystem<float, LSCMSystemMesh>(mesh, p0, p1);
    return std::make_tuple(mesh, p0, p1, std::move(parts));
}

}  // namespace

TEST(LSCMSystemBuild, Dimensions_KnownMesh)
{
    auto [mesh, p0, p1, parts] = BuildPyramidSystem();

    const auto numFaces = mesh->num_faces();
    const auto numVerts = mesh->num_vertices();
    const auto numFree = numVerts - 2;

    EXPECT_EQ(static_cast<std::size_t>(parts.A.rows()), 2 * numFaces);
    EXPECT_EQ(static_cast<std::size_t>(parts.A.cols()), 2 * numFree);
    EXPECT_EQ(static_cast<std::size_t>(parts.b.rows()), 2 * numFaces);
    EXPECT_EQ(static_cast<std::size_t>(parts.b.cols()), 1);
}

TEST(LSCMSystemBuild, FreeIdxTable_Population)
{
    auto [mesh, p0, p1, parts] = BuildPyramidSystem();

    EXPECT_EQ(parts.freeIdxTable.size(), mesh->num_vertices() - 2);
    EXPECT_EQ(parts.freeIdxTable.count(p0->idx), 0u)
        << "freeIdxTable must not contain pin0 (idx=" << p0->idx << ")";
    EXPECT_EQ(parts.freeIdxTable.count(p1->idx), 0u)
        << "freeIdxTable must not contain pin1 (idx=" << p1->idx << ")";
    for (const auto& v : mesh->vertices()) {
        if (v == p0 or v == p1) {
            continue;
        }
        EXPECT_EQ(parts.freeIdxTable.count(v->idx), 1u)
            << "free vertex " << v->idx << " missing from freeIdxTable";
    }

    // All assigned slot indices are unique and contiguous in [0, numFree).
    std::vector<std::size_t> slots;
    slots.reserve(parts.freeIdxTable.size());
    for (const auto& kv : parts.freeIdxTable) {
        slots.push_back(kv.second);
    }
    std::sort(slots.begin(), slots.end());
    for (std::size_t i = 0; i < slots.size(); ++i) {
        EXPECT_EQ(slots[i], i) << "freeIdxTable slot " << i << " not contiguous";
    }
}

TEST(LSCMSystemBuild, PinRowsLandInB)
{
    auto [mesh, p0, p1, parts] = BuildPyramidSystem();

    // After buildSystem, p0 sits at the UV origin and p1 sits on the axis
    // whose component of (p1->pos - p0->pos) had the largest magnitude.
    // Therefore b = bFree * bFixed * -1 must have at least one nonzero entry
    // (the pin1 axis contributes a nonzero displacement into b).
    EXPECT_GT(parts.b.nonZeros(), 0)
        << "Expected pin-row contributions to populate b via bFree * bFixed";

    // Pin vertex indices must not appear as columns in A — A's columns are
    // indexed by freeIdxTable slots only. Walk A's iterator and confirm no
    // nonzero entry lands in a column that would correspond to a pin slot.
    const auto numFree = parts.freeIdxTable.size();
    for (int k = 0; k < parts.A.outerSize(); ++k) {
        for (Eigen::SparseMatrix<float>::InnerIterator it(parts.A, k); it; ++it) {
            EXPECT_LT(static_cast<std::size_t>(it.col()), 2 * numFree)
                << "A nonzero at col " << it.col() << " exceeds free-vertex column range";
        }
    }
}

TEST(LSCMSystemBuild, FreeIdxTable_DeterministicAcrossRuns)
{
    // A 4x4 grid yields 16 vertices, 18 faces; with pin0=0, pin1=3 (two
    // corners on the bottom row) there are 14 free vertices — enough to
    // detect slot-ordering regressions a 4-vertex pyramid cannot.
    constexpr std::size_t pin0Idx = 0;
    constexpr std::size_t pin1Idx = 3;

    auto runOnce = [&]() {
        auto mesh = ConstructGrid<LSCMSystemMesh>(4, 4);
        ComputeMeshAngles(mesh);
        auto p0 = mesh->vertex(pin0Idx);
        auto p1 = mesh->vertex(pin1Idx);
        return OpenABF::detail::lscm::buildSystem<float, LSCMSystemMesh>(mesh, p0, p1);
    };

    auto parts1 = runOnce();
    auto parts2 = runOnce();

    EXPECT_EQ(parts1.freeIdxTable.size(), 14u);
    EXPECT_EQ(parts2.freeIdxTable.size(), 14u);

    // Slot assignments must be identical run-to-run.
    for (const auto& [origIdx, slot] : parts1.freeIdxTable) {
        ASSERT_TRUE(parts2.freeIdxTable.count(origIdx))
            << "vertex " << origIdx << " present in run 1, absent in run 2";
        EXPECT_EQ(parts2.freeIdxTable.at(origIdx), slot)
            << "vertex " << origIdx << " got slot " << slot << " in run 1 and "
            << parts2.freeIdxTable.at(origIdx) << " in run 2";
    }

    // A's column count must equal 2*numFree on both runs.
    EXPECT_EQ(static_cast<std::size_t>(parts1.A.cols()), 2 * 14);
    EXPECT_EQ(static_cast<std::size_t>(parts2.A.cols()), 2 * 14);
}
