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

#include "OpenABF/HalfEdgeMesh.hpp"

namespace OpenABF
{

/**
 * @brief Run an angle optimizer and parameterizer on every connected component
 * of @p mesh, then write the resulting UV coordinates back onto the source.
 *
 * For each connected component:
 *   1. Extract an independent mesh via `mesh->extract_connected_components()`.
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
    auto components = mesh->extract_connected_components();
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
