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
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

#include "OpenABF/Vec.hpp"

namespace OpenABF
{

namespace detail
{
/** @brief Deduce the dimensionality of a Vec type */
template <typename>
struct VecDimensions;

template <typename U, std::size_t N>
struct VecDimensions<Vec<U, N>> {
    static constexpr std::size_t value = N;
};
}  // namespace detail

/**
 * @brief Options controlling PackCharts behavior
 *
 * @tparam T Floating-point scalar type
 */
template <typename T>
struct PackOptions {
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
     * yields a roughly square atlas.
     */
    std::optional<T> target_width{};

    /**
     * @brief Gutter added around every chart, in chart/absolute units
     *
     * Applied on all four sides of each chart, including against the atlas
     * boundary, so perimeter charts are inset from the returned extent by
     * `padding` as well -- not merely separated from their neighbors.
     * Defaults to `0` (charts laid out flush). `padding` is in absolute chart
     * units and is applied before any `normalize` scaling.
     */
    T padding{T(0)};
};

/**
 * @brief The bounding box of the packed atlas
 *
 * @tparam T Floating-point scalar type
 */
template <typename T>
struct PackResult {
    /** @brief Lower corner of the packed atlas */
    Vec<T, 2> min;
    /** @brief Upper corner of the packed atlas */
    Vec<T, 2> max;
};

/**
 * @brief Pack a set of parameterized charts into a shared coordinate frame
 *
 * Lays out a list of already-parameterized charts (2D meshes whose vertex
 * `pos` holds `{u, v, ...}`) into a single shared frame using shelf packing,
 * so that no two charts' bounding boxes overlap. Operates purely on geometry:
 * each chart's vertex positions are translated (and, when `normalize` is set,
 * uniformly scaled) **in place**. Topology, vertex indices, and face indices
 * are untouched, so any `ExtractedComponent` back-maps a caller holds remain
 * valid after packing.
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
 * boxes, one to apply the transform). Memory overhead is `O(n)`.
 *
 * @tparam MeshType A HalfEdgeMesh specialization
 * @param charts Charts to pack; each chart's vertex positions are modified
 * @param opts Packing options
 * @return The bounding box of the packed atlas
 *
 * @throws std::invalid_argument If a chart pointer is null or a chart has no
 *         vertices.
 */
template <typename MeshType>
auto PackCharts(std::vector<typename MeshType::Pointer>& charts,
                const PackOptions<typename MeshType::type>& opts =
                    PackOptions<typename MeshType::type>{}) -> PackResult<typename MeshType::type>
{
    using T = typename MeshType::type;
    static_assert(
        detail::VecDimensions<decltype(std::declval<typename MeshType::Vertex>().pos)>::value >= 2,
        "PackCharts requires mesh vertices with at least 2 position dimensions");

    PackResult<T> result{Vec<T, 2>{T(0), T(0)}, Vec<T, 2>{T(0), T(0)}};
    if (charts.empty()) {
        return result;
    }
    const std::size_t n = charts.size();

    // Validate inputs and compute each chart's 2D bounding box (origin + size).
    std::vector<T> minX(n);
    std::vector<T> minY(n);
    std::vector<T> width(n);
    std::vector<T> height(n);
    for (std::size_t i = 0; i < n; ++i) {
        const auto& chart = charts[i];
        if (not chart or chart->num_vertices() == 0) {
            throw std::invalid_argument("PackCharts: chart is null or has no vertices");
        }
        auto mnX = std::numeric_limits<T>::max();
        auto mnY = std::numeric_limits<T>::max();
        auto mxX = std::numeric_limits<T>::lowest();
        auto mxY = std::numeric_limits<T>::lowest();
        for (const auto& v : chart->vertices()) {
            mnX = std::min(mnX, v->pos[0]);
            mnY = std::min(mnY, v->pos[1]);
            mxX = std::max(mxX, v->pos[0]);
            mxY = std::max(mxY, v->pos[1]);
        }
        minX[i] = mnX;
        minY[i] = mnY;
        width[i] = mxX - mnX;
        height[i] = mxY - mnY;
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
    auto atlasMaxX = T(0);
    auto atlasMaxY = T(0);
    for (const auto i : order) {
        if (cursorX > pad and cursorX + width[i] > targetWidth) {
            cursorX = pad;
            cursorY += shelfHeight + pad;
            shelfHeight = T(0);
        }
        offsetX[i] = cursorX - minX[i];
        offsetY[i] = cursorY - minY[i];
        atlasMaxX = std::max(atlasMaxX, cursorX + width[i]);
        atlasMaxY = std::max(atlasMaxY, cursorY + height[i]);
        cursorX += width[i] + pad;
        shelfHeight = std::max(shelfHeight, height[i]);
    }

    // The atlas extent includes the perimeter gutter: charts are inset by
    // `pad` from the lower corner (the cursor starts at `pad`), so add `pad`
    // to the far edges too. Every chart then has >= `pad` of empty space on
    // all four sides, including against the atlas boundary.
    const T atlasW = atlasMaxX + pad;
    const T atlasH = atlasMaxY + pad;

    // Optional normalization: a single global uniform scale that fits the
    // padded atlas into [0,1]^2, preserving relative chart sizes.
    auto scale = T(1);
    if (opts.normalize) {
        const auto extentMax = std::max(atlasW, atlasH);
        if (extentMax > T(0)) {
            scale = T(1) / extentMax;
        }
    }

    // Apply translation (+ optional scale about the origin) in place.
    for (std::size_t i = 0; i < n; ++i) {
        for (const auto& v : charts[i]->vertices()) {
            v->pos[0] = (v->pos[0] + offsetX[i]) * scale;
            v->pos[1] = (v->pos[1] + offsetY[i]) * scale;
        }
    }

    result.max = Vec<T, 2>{atlasW * scale, atlasH * scale};
    return result;
}

}  // namespace OpenABF
