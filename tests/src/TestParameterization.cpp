#include <cmath>
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