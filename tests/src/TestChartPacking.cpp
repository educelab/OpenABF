#include <algorithm>
#include <limits>
#include <set>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "OpenABF/OpenABF.hpp"

using Mesh = OpenABF::HalfEdgeMesh<float>;
using OpenABF::PackCharts;
using OpenABF::PackOptions;
using OpenABF::PackResult;

namespace
{

/** @brief Build a w x h rectangle chart (two triangles) at the origin */
auto MakeRectChart(float w, float h) -> Mesh::Pointer
{
    auto m = Mesh::New();
    m->insert_vertices({{0.f, 0.f, 0.f}, {w, 0.f, 0.f}, {w, h, 0.f}, {0.f, h, 0.f}});
    m->insert_faces({{0, 1, 2}, {0, 2, 3}});
    return m;
}

struct BBox {
    float minx{std::numeric_limits<float>::max()};
    float miny{std::numeric_limits<float>::max()};
    float maxx{std::numeric_limits<float>::lowest()};
    float maxy{std::numeric_limits<float>::lowest()};
    [[nodiscard]] auto width() const -> float { return maxx - minx; }
    [[nodiscard]] auto height() const -> float { return maxy - miny; }
};

template <typename MeshPtr>
auto ChartBBox(const MeshPtr& m) -> BBox
{
    BBox b;
    for (const auto& v : m->vertices()) {
        b.minx = std::min(b.minx, v->pos[0]);
        b.miny = std::min(b.miny, v->pos[1]);
        b.maxx = std::max(b.maxx, v->pos[0]);
        b.maxy = std::max(b.maxy, v->pos[1]);
    }
    return b;
}

/** @brief True if two boxes overlap by more than eps (touching is allowed) */
auto Overlaps(const BBox& a, const BBox& b, float eps) -> bool
{
    return (a.minx + eps < b.maxx) && (b.minx + eps < a.maxx) && (a.miny + eps < b.maxy) &&
           (b.miny + eps < a.maxy);
}

}  // namespace

// --- Edge cases ------------------------------------------------------------

TEST(ChartPacking, EmptyListReturnsEmptyExtent)
{
    std::vector<Mesh::Pointer> charts;
    auto extent = PackCharts<Mesh>(charts);
    EXPECT_FLOAT_EQ(extent.min[0], 0.f);
    EXPECT_FLOAT_EQ(extent.min[1], 0.f);
    EXPECT_FLOAT_EQ(extent.max[0], 0.f);
    EXPECT_FLOAT_EQ(extent.max[1], 0.f);
}

TEST(ChartPacking, NullChartThrows)
{
    std::vector<Mesh::Pointer> charts{MakeRectChart(1.f, 1.f), nullptr};
    EXPECT_THROW(PackCharts<Mesh>(charts), std::invalid_argument);
}

TEST(ChartPacking, EmptyChartThrows)
{
    std::vector<Mesh::Pointer> charts{MakeRectChart(1.f, 1.f), Mesh::New()};
    EXPECT_THROW(PackCharts<Mesh>(charts), std::invalid_argument);
}

TEST(ChartPacking, ZeroAreaChartIsPlacedWithoutThrowing)
{
    // A flat (collinear) chart has zero area but is a valid mesh: distinct
    // vertices, non-zero edges. PackCharts must place it without crashing.
    auto flat = Mesh::New();
    flat->insert_vertices({{0.f, 0.f, 0.f}, {1.f, 0.f, 0.f}, {2.f, 0.f, 0.f}});
    flat->insert_faces({{0, 1, 2}});
    std::vector<Mesh::Pointer> charts{MakeRectChart(1.f, 1.f), flat};
    EXPECT_NO_THROW(PackCharts<Mesh>(charts));
}

// --- Absolute (default) mode ----------------------------------------------

TEST(ChartPacking, SingleChartMapsMinToOrigin)
{
    std::vector<Mesh::Pointer> charts{MakeRectChart(3.f, 2.f)};
    auto extent = PackCharts<Mesh>(charts);
    auto b = ChartBBox(charts[0]);
    EXPECT_FLOAT_EQ(b.minx, 0.f);
    EXPECT_FLOAT_EQ(b.miny, 0.f);
    EXPECT_FLOAT_EQ(b.width(), 3.f);
    EXPECT_FLOAT_EQ(b.height(), 2.f);
    EXPECT_FLOAT_EQ(extent.max[0], 3.f);
    EXPECT_FLOAT_EQ(extent.max[1], 2.f);
}

TEST(ChartPacking, AbsoluteModePreservesChartSizes)
{
    std::vector<float> ws{1.f, 2.f, 0.5f, 3.f};
    std::vector<float> hs{1.f, 1.5f, 2.f, 0.7f};
    std::vector<Mesh::Pointer> charts;
    for (std::size_t i = 0; i < ws.size(); ++i) {
        charts.push_back(MakeRectChart(ws[i], hs[i]));
    }
    PackCharts<Mesh>(charts);
    for (std::size_t i = 0; i < charts.size(); ++i) {
        auto b = ChartBBox(charts[i]);
        EXPECT_FLOAT_EQ(b.width(), ws[i]) << "chart " << i;
        EXPECT_FLOAT_EQ(b.height(), hs[i]) << "chart " << i;
    }
}

TEST(ChartPacking, ChartsDoNotOverlap)
{
    std::vector<Mesh::Pointer> charts;
    for (int i = 0; i < 6; ++i) {
        charts.push_back(MakeRectChart(1.f + 0.3f * static_cast<float>(i), 1.f));
    }
    PackCharts<Mesh>(charts);
    std::vector<BBox> boxes;
    for (const auto& c : charts) {
        boxes.push_back(ChartBBox(c));
    }
    for (std::size_t i = 0; i < boxes.size(); ++i) {
        for (std::size_t j = i + 1; j < boxes.size(); ++j) {
            EXPECT_FALSE(Overlaps(boxes[i], boxes[j], 1e-4f)) << "charts " << i << " and " << j;
        }
    }
}

TEST(ChartPacking, ExtentBoundsAllCharts)
{
    std::vector<Mesh::Pointer> charts{MakeRectChart(1.f, 1.f), MakeRectChart(2.f, 0.5f),
                                      MakeRectChart(0.5f, 3.f)};
    auto extent = PackCharts<Mesh>(charts);
    for (const auto& c : charts) {
        for (const auto& v : c->vertices()) {
            EXPECT_GE(v->pos[0], extent.min[0] - 1e-4f);
            EXPECT_GE(v->pos[1], extent.min[1] - 1e-4f);
            EXPECT_LE(v->pos[0], extent.max[0] + 1e-4f);
            EXPECT_LE(v->pos[1], extent.max[1] + 1e-4f);
        }
    }
}

TEST(ChartPacking, PaddingSeparatesChartsInSingleRow)
{
    // Force a single row with a large target width; padding should appear as a
    // gap between the two charts, so the atlas width is w0 + padding + w1.
    std::vector<Mesh::Pointer> charts{MakeRectChart(1.f, 1.f), MakeRectChart(1.f, 1.f)};
    PackOptions<float> opts;
    opts.target_width = 1000.f;
    opts.padding = 0.5f;
    auto extent = PackCharts<Mesh>(charts, opts);
    EXPECT_NEAR(extent.max[0] - extent.min[0], 2.5f, 1e-4f);

    auto b0 = ChartBBox(charts[0]);
    auto b1 = ChartBBox(charts[1]);
    auto gap = std::max(b1.minx - b0.maxx, b0.minx - b1.maxx);
    EXPECT_GE(gap, 0.5f - 1e-4f);
}

// --- Normalize mode --------------------------------------------------------

TEST(ChartPacking, NormalizeFitsUnitSquare)
{
    std::vector<Mesh::Pointer> charts{MakeRectChart(2.f, 1.f), MakeRectChart(1.f, 3.f),
                                      MakeRectChart(0.5f, 0.5f)};
    PackOptions<float> opts;
    opts.normalize = true;
    auto extent = PackCharts<Mesh>(charts, opts);

    float maxCoord = 0.f;
    for (const auto& c : charts) {
        for (const auto& v : c->vertices()) {
            EXPECT_GE(v->pos[0], -1e-4f);
            EXPECT_GE(v->pos[1], -1e-4f);
            EXPECT_LE(v->pos[0], 1.f + 1e-4f);
            EXPECT_LE(v->pos[1], 1.f + 1e-4f);
            maxCoord = std::max({maxCoord, v->pos[0], v->pos[1]});
        }
    }
    // The atlas should fill at least one axis of the unit square.
    EXPECT_NEAR(maxCoord, 1.f, 1e-3f);
    EXPECT_LE(extent.max[0], 1.f + 1e-4f);
    EXPECT_LE(extent.max[1], 1.f + 1e-4f);
}

TEST(ChartPacking, NormalizePreservesRelativeChartSizes)
{
    // One chart twice the linear size of the other: ratio must survive a single
    // global uniform scale.
    std::vector<Mesh::Pointer> charts{MakeRectChart(1.f, 1.f), MakeRectChart(2.f, 2.f)};
    PackOptions<float> opts;
    opts.normalize = true;
    PackCharts<Mesh>(charts, opts);
    auto b0 = ChartBBox(charts[0]);
    auto b1 = ChartBBox(charts[1]);
    EXPECT_NEAR(b1.width() / b0.width(), 2.f, 1e-3f);
    EXPECT_NEAR(b1.height() / b0.height(), 2.f, 1e-3f);
}

// --- End-to-end pipeline + per-wedge recovery ------------------------------

TEST(ChartPacking, EndToEndTearExtractFlattenPackAndWedgeRecovery)
{
    using ABF = OpenABF::ABFPlusPlus<float>;
    using LSCM = OpenABF::AngleBasedLSCM<float, ABF::Mesh>;

    // 3x3 grid (matches the MultiChartFlatten example).
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

    // Tear into two charts and extract.
    mesh->split_path({1, 4, 7});
    auto ccs = mesh->extract_connected_components();
    ASSERT_EQ(ccs.size(), 2u);

    // Flatten each chart and collect the meshes for packing.
    std::vector<ABF::Mesh::Pointer> chartMeshes;
    for (auto& cc : ccs) {
        std::size_t iters{0};
        float grad{OpenABF::INF<float>};
        ABF::Compute(cc.mesh, iters, grad);
        LSCM::Compute(cc.mesh);
        chartMeshes.push_back(cc.mesh);
    }

    auto extent = PackCharts<ABF::Mesh>(chartMeshes);

    // Charts must not overlap after packing.
    std::vector<BBox> boxes;
    for (const auto& c : chartMeshes) {
        boxes.push_back(ChartBBox(c));
    }
    EXPECT_FALSE(Overlaps(boxes[0], boxes[1], 1e-4f));

    // Per-wedge recovery via the back-maps, keyed on vertex identity. Every
    // (original face, original vertex) wedge must be unique, proving the
    // documented recipe yields a well-formed per-wedge map.
    std::set<std::pair<std::size_t, std::size_t>> wedges;
    std::size_t corners = 0;
    for (auto& cc : ccs) {
        for (const auto& face : cc.mesh->faces()) {
            for (const auto& edge : *face) {
                auto origFace = cc.face_map[face->idx];
                auto origVert = cc.vertex_map[edge->vertex->idx];
                wedges.emplace(origFace, origVert);
                ++corners;
            }
        }
    }
    EXPECT_EQ(corners, 8u * 3u);  // 8 faces, 3 corners each
    EXPECT_EQ(wedges.size(), corners);
}
