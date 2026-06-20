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

#include <cstddef>
#include <optional>
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

    /** @brief Gutter added around each chart, in chart/absolute units */
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

    // STUB (TDD red phase): real layout is implemented in Phase 3.
    (void)charts;
    (void)opts;
    return PackResult<T>{Vec<T, 2>{T(0), T(0)}, Vec<T, 2>{T(0), T(0)}};
}

}  // namespace OpenABF
