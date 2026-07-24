/*
OpenABF
https://gitlab.com/educelab/OpenABF

Copyright 2025 EduceLab

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace OpenABF
{

namespace detail
{
/**
 * @brief Rotate a chart within its UV plane so its axis-aligned bounding box
 * has minimum area, with the larger extent vertical
 *
 * The minimum-area enclosing rectangle of a planar point set always has one
 * edge collinear with an edge of the set's convex hull, so it suffices to test
 * the orientation induced by each hull edge. The chart is then stood on its
 * long axis (larger extent vertical) so it aligns with PackCharts's
 * tallest-first shelf strategy. Vertex positions are rotated in place about the
 * origin; only the first two components are touched. Rotation preserves
 * topology and vertex identity, so any back-maps remain valid.
 *
 * @tparam MeshType A HalfEdgeMesh specialization
 */
template <typename MeshType>
void MinimizeChartBoundingBox(const typename MeshType::Pointer& chart)
{
    using T = typename MeshType::type;
    using Point = std::array<T, 2>;

    // Gather the 2D point set.
    std::vector<Point> pts;
    pts.reserve(chart->num_vertices());
    for (const auto& v : chart->vertices()) {
        pts.push_back({v->pos[0], v->pos[1]});
    }

    // Convex hull via Andrew's monotone chain. Fewer than three unique points
    // means a point or a segment, for which no rotation reduces the area.
    std::sort(pts.begin(), pts.end());
    pts.erase(std::unique(pts.begin(), pts.end()), pts.end());
    const std::size_t m = pts.size();
    if (m < 3) {
        return;
    }
    auto crossZ = [](const Point& o, const Point& a, const Point& b) -> T {
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
    };
    std::vector<Point> hull(2 * m);
    std::size_t k = 0;
    for (std::size_t i = 0; i < m; ++i) {
        while (k >= 2 and crossZ(hull[k - 2], hull[k - 1], pts[i]) <= T(0)) {
            --k;
        }
        hull[k++] = pts[i];
    }
    for (std::size_t i = m - 1, t = k + 1; i > 0; --i) {
        while (k >= t and crossZ(hull[k - 2], hull[k - 1], pts[i - 1]) <= T(0)) {
            --k;
        }
        hull[k++] = pts[i - 1];
    }
    hull.resize(k - 1);  // drop the duplicated start point
    const std::size_t h = hull.size();
    if (h < 3) {
        return;
    }

    // Bounding-box dimensions {width, height} of the hull after rotating every
    // point by R(-theta), where (c, s) = (cos theta, sin theta).
    auto boxFor = [&](T c, T s) -> Point {
        auto mnX = std::numeric_limits<T>::max();
        auto mnY = std::numeric_limits<T>::max();
        auto mxX = std::numeric_limits<T>::lowest();
        auto mxY = std::numeric_limits<T>::lowest();
        for (const auto& p : hull) {
            const auto rx = c * p[0] + s * p[1];
            const auto ry = -s * p[0] + c * p[1];
            mnX = std::min(mnX, rx);
            mnY = std::min(mnY, ry);
            mxX = std::max(mxX, rx);
            mxY = std::max(mxY, ry);
        }
        return {mxX - mnX, mxY - mnY};
    };

    // Seed with the current (unrotated) box so we only rotate on a strict
    // area improvement.
    auto bestCos = T(1);
    auto bestSin = T(0);
    auto bestBox = boxFor(T(1), T(0));
    auto bestArea = bestBox[0] * bestBox[1];
    for (std::size_t i = 0; i < h; ++i) {
        const auto& p0 = hull[i];
        const auto& p1 = hull[(i + 1) % h];
        const auto ex = p1[0] - p0[0];
        const auto ey = p1[1] - p0[1];
        const auto len = std::sqrt(ex * ex + ey * ey);
        if (len <= T(0)) {
            continue;
        }
        // Align this hull edge with the x-axis (theta = atan2(ey, ex)).
        const auto c = ex / len;
        const auto s = ey / len;
        const auto box = boxFor(c, s);
        const auto area = box[0] * box[1];
        if (area < bestArea) {
            bestArea = area;
            bestCos = c;
            bestSin = s;
            bestBox = box;
        }
    }

    // Stand the chart on its long axis: the packer sorts tallest-first and
    // fills horizontal shelves, so the larger extent should be vertical. If the
    // min-area box is wider than tall, compose an extra 90-degree rotation
    // (R90 * R(-theta), where R90 maps (x, y) -> (-y, x)).
    if (bestBox[0] > bestBox[1]) {
        const auto c = bestCos;
        const auto s = bestSin;
        bestCos = s;
        bestSin = -c;
    }

    // Apply the chosen rotation to every vertex in place (skip the identity).
    if (bestCos != T(1) or bestSin != T(0)) {
        for (const auto& v : chart->vertices()) {
            const auto x = v->pos[0];
            const auto y = v->pos[1];
            v->pos[0] = bestCos * x + bestSin * y;
            v->pos[1] = -bestSin * x + bestCos * y;
        }
    }
}
}  // namespace detail

/**
 * @brief Options controlling PackCharts behavior
 *
 * @tparam T Floating-point scalar type
 */
template <typename T>
struct PackOptions {
    /**
     * @brief Rotate each chart in-plane to minimize its bounding-box area
     *
     * When `true` (default), each chart is rotated within its UV plane before
     * layout so its axis-aligned bounding box has minimum area, then stood on
     * its long axis (larger extent vertical) to match the tallest-first shelf
     * strategy. Shelf packing works on axis-aligned boxes, so tightening and
     * consistently orienting each box lets charts nest more densely. The
     * rotation is applied in place and preserves topology and vertex identity,
     * so any back-maps a caller holds remain valid.
     */
    bool minimize_bounding_box{true};

    /**
     * @brief Fit the packed atlas into the unit square `[0,1]^2`
     *
     * When `false` (default), charts keep their absolute (physical) scale and
     * are only translated. When `true`, a single global uniform scale is
     * applied after layout so the whole atlas fits in `[0,1]^2`. A single
     * global factor preserves relative chart sizes and cross-chart texel
     * density; only the absolute units change.
     */
    bool normalize{false};

    /**
     * @brief Target shelf width; overrides the `sqrt(total area)` heuristic
     *
     * Charts wrap to a new shelf when a row would exceed this width. If unset,
     * the width defaults to `sqrt(sum of per-chart bounding-box areas)`, which
     * yields a roughly square atlas. Must be positive when set.
     */
    std::optional<T> target_width{};

    /**
     * @brief Gutter added around every chart, in chart/absolute units
     *
     * Applied on all four sides of each chart, including against the atlas
     * boundary, so perimeter charts are inset from the returned extent by
     * `padding` as well -- not merely separated from their neighbors.
     * Defaults to `0` (charts laid out flush). `padding` is in absolute chart
     * units and is applied before any `normalize` scaling. Must be
     * non-negative.
     */
    T padding{T(0)};
};

/**
 * @brief The bounding box of the packed atlas
 *
 * `min`/`max` use the same vector type as the input meshes' vertex positions;
 * only the first two (`u`, `v`) components are meaningful.
 *
 * @tparam VecType The vertex position vector type of the packed meshes
 */
template <typename VecType>
struct PackResult {
    /** @brief Lower corner of the packed atlas */
    VecType min;
    /** @brief Upper corner of the packed atlas */
    VecType max;
};

/**
 * @brief Pack a set of parameterized charts into a shared coordinate frame
 *
 * Lays out a list of already-parameterized charts (2D meshes whose vertex
 * `pos` holds `{u, v, ...}`) into a single shared frame using shelf packing,
 * so that no two charts' bounding boxes overlap. Operates purely on geometry:
 * each chart's vertex positions are rotated (when `minimize_bounding_box` is
 * set), translated, and (when `normalize` is set) uniformly scaled **in
 * place**. Topology, vertex indices, and face indices are untouched, so any
 * `ExtractedComponent` back-maps a caller holds remain valid after packing.
 *
 * @par Scaling
 * By default charts keep their absolute scale and are only translated; the
 * returned extent is meaningful in physical units. With `opts.normalize`, one
 * global uniform scale maps the packed atlas into `[0,1]^2`.
 *
 * @par Building a per-wedge UV map
 * This function does not own a UV-map type. To build a per-corner ("wedge")
 * UV map from the packed charts, key each wedge by **vertex identity**, not by
 * corner position: a face's corner order is not stable (the half-edge mesh may
 * reverse a mis-wound face at insertion time, and that permutation is not
 * recorded). Given an `ExtractedComponent` `ec` for a chart, the robust key is:
 * @code
 * for (const auto& face : chart->faces()) {
 *     for (const auto& edge : *face) {            // walk the face's corners
 *         auto origFace = ec.face_map[face->idx];
 *         auto origVert = ec.vertex_map[edge->vertex->idx];
 *         auto uv       = edge->vertex->pos;       // packed UV
 *         // wedge (origFace, origVert) -> uv
 *     }
 * }
 * @endcode
 *
 * @par Complexity
 * `O(n log n)` in the number of charts `n` (dominated by the height sort) plus
 * `O(V)` in the total vertex count `V` (two passes: one to measure bounding
 * boxes, one to apply the transform). With `minimize_bounding_box`, each chart
 * additionally costs an `O(v log v)` convex hull over its `v` vertices plus an
 * `O(h^2)` orientation search, which measures the hull's `h` points once per
 * hull edge. Memory overhead is `O(n)`.
 *
 * @tparam MeshType A HalfEdgeMesh specialization
 * @param charts Charts to pack; each chart's vertex positions are modified
 * @param opts Packing options
 * @return The bounding box of the packed atlas
 *
 * @throws std::invalid_argument If a chart pointer is null, a chart has no
 *         vertices, `opts.padding` is negative, or `opts.target_width` is set
 *         to a non-positive value. All inputs are validated before any chart is
 *         modified, so a throw leaves every chart untouched.
 */
template <typename MeshType>
auto PackCharts(
    std::vector<std::shared_ptr<MeshType>>& charts,
    const PackOptions<typename MeshType::type>& opts = PackOptions<typename MeshType::type>{})
    -> PackResult<typename MeshType::PositionType>
{
    using T = typename MeshType::type;
    using VecType = typename MeshType::PositionType;
    static_assert(VecType::Dimensions >= 2,
                  "PackCharts requires mesh vertices with at least 2 position dimensions");

    PackResult<VecType> result{};
    if (charts.empty()) {
        return result;
    }
    const std::size_t n = charts.size();

    // Validate everything up front. This function mutates caller-owned meshes in
    // place, so throwing partway through would leave some charts rotated and
    // others not, with no way for the caller to tell which.
    if (opts.padding < T(0)) {
        throw std::invalid_argument("PackCharts: padding must be non-negative");
    }
    if (opts.target_width and *opts.target_width <= T(0)) {
        throw std::invalid_argument("PackCharts: target_width must be positive");
    }
    for (const auto& chart : charts) {
        if (not chart or chart->num_vertices() == 0) {
            throw std::invalid_argument("PackCharts: chart is null or has no vertices");
        }
    }

    // Compute each chart's 2D bounding box (origin + size).
    std::vector<T> minX(n);
    std::vector<T> minY(n);
    std::vector<T> width(n);
    std::vector<T> height(n);
    for (std::size_t i = 0; i < n; ++i) {
        const auto& chart = charts[i];
        // Tighten the chart's bounding box by rotating it in-plane first.
        if (opts.minimize_bounding_box) {
            detail::MinimizeChartBoundingBox<MeshType>(chart);
        }
        const auto& verts = chart->vertices();
        const auto cmpX = [](const auto& a, const auto& b) { return a->pos[0] < b->pos[0]; };
        const auto cmpY = [](const auto& a, const auto& b) { return a->pos[1] < b->pos[1]; };
        const auto [xlo, xhi] = std::minmax_element(verts.begin(), verts.end(), cmpX);
        const auto [ylo, yhi] = std::minmax_element(verts.begin(), verts.end(), cmpY);
        minX[i] = (*xlo)->pos[0];
        minY[i] = (*ylo)->pos[1];
        width[i] = (*xhi)->pos[0] - (*xlo)->pos[0];
        height[i] = (*yhi)->pos[1] - (*ylo)->pos[1];
    }

    // Place taller charts first so shelves pack tightly.
    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::sort(order.begin(), order.end(),
              [&](std::size_t a, std::size_t b) { return height[a] > height[b]; });

    // Target shelf width: caller override, else sqrt(sum of padded chart areas)
    // for a roughly square atlas.
    const T pad = opts.padding;
    T targetWidth;
    if (opts.target_width) {
        targetWidth = *opts.target_width;
    } else {
        auto areaSum = T(0);
        for (std::size_t i = 0; i < n; ++i) {
            areaSum += (width[i] + pad) * (height[i] + pad);
        }
        targetWidth = std::sqrt(areaSum);
    }

    // Shelf layout. Charts are placed left-to-right; a row wraps to a new shelf
    // once it would exceed targetWidth (a chart wider than targetWidth still
    // gets placed alone at the start of a shelf). The cursor starts at `pad`
    // and wraps back to `pad`, so every chart is inset by at least `pad` from
    // the atlas's lower corner; the atlas's lower corner itself stays at the
    // origin.
    std::vector<T> offsetX(n);
    std::vector<T> offsetY(n);
    auto cursorX = pad;
    auto cursorY = pad;
    auto shelfHeight = T(0);
    auto placedMinX = std::numeric_limits<T>::max();
    auto placedMinY = std::numeric_limits<T>::max();
    auto placedMaxX = std::numeric_limits<T>::lowest();
    auto placedMaxY = std::numeric_limits<T>::lowest();
    for (const auto i : order) {
        if (cursorX > pad and cursorX + width[i] > targetWidth) {
            cursorX = pad;
            cursorY += shelfHeight + pad;
            shelfHeight = T(0);
        }
        offsetX[i] = cursorX - minX[i];
        offsetY[i] = cursorY - minY[i];
        placedMinX = std::min(placedMinX, cursorX);
        placedMinY = std::min(placedMinY, cursorY);
        placedMaxX = std::max(placedMaxX, cursorX + width[i]);
        placedMaxY = std::max(placedMaxY, cursorY + height[i]);
        cursorX += width[i] + pad;
        shelfHeight = std::max(shelfHeight, height[i]);
    }

    // The atlas extent is measured from where the charts actually landed, then
    // grown by the perimeter gutter, so every chart has >= `pad` of empty space
    // on all four sides including against the atlas boundary. Measuring rather
    // than assuming the lower corner keeps the returned extent honest if the
    // layout strategy above ever changes.
    const T atlasMinX = placedMinX - pad;
    const T atlasMinY = placedMinY - pad;
    const T atlasW = (placedMaxX + pad) - atlasMinX;
    const T atlasH = (placedMaxY + pad) - atlasMinY;

    // Optional normalization: a single global uniform scale that fits the
    // padded atlas into [0,1]^2, preserving relative chart sizes.
    auto scale = T(1);
    if (opts.normalize) {
        const auto extentMax = std::max(atlasW, atlasH);
        if (extentMax > T(0)) {
            scale = T(1) / extentMax;
        }
    }

    // Apply translation (+ optional scale about the origin) in place. Shifting
    // by the measured atlas corner puts that corner on the origin by
    // construction rather than by assumption; with the shelf layout above the
    // shift is already zero.
    for (std::size_t i = 0; i < n; ++i) {
        for (const auto& v : charts[i]->vertices()) {
            v->pos[0] = (v->pos[0] + offsetX[i] - atlasMinX) * scale;
            v->pos[1] = (v->pos[1] + offsetY[i] - atlasMinY) * scale;
        }
    }

    // The atlas lower corner now sits on the origin, so result.min keeps its
    // value-initialized zero; only u/v of max are set.
    result.max[0] = atlasW * scale;
    result.max[1] = atlasH * scale;
    return result;
}

}  // namespace OpenABF
