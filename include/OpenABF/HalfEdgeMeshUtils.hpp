/*
OpenABF
https://gitlab.com/educelab/OpenABF

Copyright 2025 EduceLab

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0
*/
#pragma once

#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#include "OpenABF/HalfEdgeMesh.hpp"

namespace OpenABF
{

/**
 * @brief Extract each connected component of @p mesh into its own independent
 * mesh.
 *
 * The returned vector contains one `(extracted, back_map)` pair per connected
 * component. The extracted mesh is a deep copy: vertices are inserted via the
 * Vertex copy ctor (so user-defined VertexTraits are preserved), faces are
 * inserted via `insert_faces` which recomputes per-face geometric angles
 * (`alpha`) from positions. Vertex indices in the extracted mesh are
 * re-densified to `0..N-1` and `back_map[extracted_idx] == original_idx`.
 *
 * A single-component mesh round-trips as a single extracted copy whose
 * back-map is `{0, 1, ..., N-1}`.
 *
 * @note Edge traits beyond the geometric `alpha` are not currently copied; if
 * your code stores per-edge state that is not derivable from vertex positions,
 * recompute it on the extracted mesh after this call.
 */
template <typename Mesh>
auto ExtractConnectedComponents(const typename Mesh::Pointer& mesh)
    -> std::vector<std::pair<typename Mesh::Pointer, std::vector<std::size_t>>>
{
    using MeshPtr = typename Mesh::Pointer;
    std::vector<std::pair<MeshPtr, std::vector<std::size_t>>> result;

    for (const auto& component : mesh->connected_components()) {
        auto sub = Mesh::New();
        std::unordered_map<std::size_t, std::size_t> remap;
        std::vector<std::size_t> backMap;
        std::vector<std::vector<std::size_t>> faceIdxs;
        faceIdxs.reserve(component.size());

        for (const auto& face : component) {
            std::vector<std::size_t> tri;
            for (const auto& edge : *face) {
                const auto& v = edge->vertex;
                auto it = remap.find(v->idx);
                if (it == remap.end()) {
                    const auto newIdx = sub->insert_vertex(*v);
                    it = remap.emplace(v->idx, newIdx).first;
                    backMap.push_back(v->idx);
                }
                tri.push_back(it->second);
            }
            faceIdxs.push_back(std::move(tri));
        }

        sub->insert_faces(faceIdxs);
        result.emplace_back(std::move(sub), std::move(backMap));
    }

    return result;
}

/**
 * @brief Run an angle optimizer and parameterizer on every connected component
 * of @p mesh, then write the resulting UV coordinates back onto the source.
 *
 * For each connected component:
 *   1. Extract an independent mesh via `ExtractConnectedComponents`.
 *   2. Optionally run `AngleOptimizer::Compute(sub)` (skipped if
 *      `AngleOptimizer` is `void`).
 *   3. Run `Parameterizer::Compute(sub)`.
 *   4. Copy `sub->vertex(i)->pos` back onto `mesh->vertex(back_map[i])->pos`.
 *
 * Use this when a mesh has been torn into multiple pieces (e.g. via
 * `split_path`) and you want to parameterize each piece with a single call.
 *
 * @tparam AngleOptimizer Type with a `static Compute(MeshPtr&)` method
 *         (e.g. `ABF`, `ABFPlusPlus`), or `void` to skip the angle pass.
 * @tparam Parameterizer Type with a `static Compute(MeshPtr&)` method
 *         (e.g. `AngleBasedLSCM`, `HierarchicalLSCM`).
 */
template <typename AngleOptimizer, typename Parameterizer, typename MeshPtr>
void ParameterizeConnectedComponents(const MeshPtr& mesh)
{
    using Mesh = typename MeshPtr::element_type;
    auto components = ExtractConnectedComponents<Mesh>(mesh);
    for (auto& [sub, backMap] : components) {
        if constexpr (not std::is_void_v<AngleOptimizer>) {
            AngleOptimizer::Compute(sub);
        }
        Parameterizer::Compute(sub);
        for (std::size_t i = 0; i < sub->num_vertices(); ++i) {
            mesh->vertex(backMap[i])->pos = sub->vertex(i)->pos;
        }
    }
}

}  // namespace OpenABF
