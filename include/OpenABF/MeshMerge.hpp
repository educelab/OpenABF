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
#include <stdexcept>
#include <utility>
#include <vector>

namespace OpenABF
{

/**
 * @brief Result of merging several meshes into one
 *
 * The inverse bookkeeping of `HalfEdgeMesh::extract_connected_components`:
 * where extract splits one mesh into N with maps back to the source, merge
 * combines N meshes into one with maps back to each source mesh. Compose
 * `face_source`/`vertex_source` with the `face_map`/`vertex_map` of the
 * components you merged to trace a merged element all the way back to the
 * pre-extraction mesh.
 *
 * @tparam MeshType A HalfEdgeMesh specialization
 */
template <typename MeshType>
struct MergedMesh {
    /** @brief The combined mesh */
    typename MeshType::Pointer mesh;
    /** @brief `vertex_source[merged_idx] = {input mesh index, source vertex index}` */
    std::vector<std::pair<std::size_t, std::size_t>> vertex_source;
    /** @brief `face_source[merged_idx] = {input mesh index, source face index}` */
    std::vector<std::pair<std::size_t, std::size_t>> face_source;
};

/**
 * @brief Merge several meshes into a single mesh
 *
 * Concatenates the vertices and faces of each input mesh into one new mesh,
 * offsetting face vertex indices per input so the inputs remain disjoint
 * components. Returns the combined mesh alongside `vertex_source`/`face_source`
 * maps that record, for every merged vertex and face, which input mesh and
 * which source index it came from — the inverse of
 * `extract_connected_components`, so the back-map chain survives the merge.
 *
 * Vertex positions and vertex traits are preserved (via the vertex copy
 * constructor). Edge and face traits are default-constructed: merge rebuilds
 * connectivity, so per-edge/per-face solver state is not carried over.
 *
 * @tparam MeshType A HalfEdgeMesh specialization
 * @param meshes Meshes to merge
 * @return The combined mesh and its provenance maps
 *
 * @throws std::invalid_argument If an input pointer is null or has no vertices.
 */
template <typename MeshType>
auto MergeMeshes(const std::vector<typename MeshType::Pointer>& meshes) -> MergedMesh<MeshType>
{
    MergedMesh<MeshType> result{MeshType::New(), {}, {}};
    auto& out = result.mesh;

    // Gather every face's (offset) vertex indices so they can be inserted in a
    // single insert_faces() call, which rebuilds the mesh boundary once at the
    // end via update_boundary(). Inserting faces one at a time with
    // insert_face() would leave the boundary stale.
    std::vector<std::vector<std::size_t>> faces;
    for (std::size_t ci = 0; ci < meshes.size(); ++ci) {
        const auto& src = meshes[ci];
        if (not src or src->num_vertices() == 0) {
            throw std::invalid_argument("MergeMeshes: input mesh is null or has no vertices");
        }
        // Vertices keep their relative order, so merged index == offset + sub
        // index; record provenance alongside each insertion.
        const auto offset = out->num_vertices();
        for (const auto& v : src->vertices()) {
            out->insert_vertex(*v);
            result.vertex_source.emplace_back(ci, v->idx);
        }
        // Re-emit each face against the offset vertex indices, preserving the
        // source face's corner order.
        for (const auto& face : src->faces()) {
            std::vector<std::size_t> idxs;
            for (const auto& edge : *face) {
                idxs.push_back(edge->vertex->idx + offset);
            }
            faces.push_back(std::move(idxs));
            result.face_source.emplace_back(ci, face->idx);
        }
    }
    out->insert_faces(faces);
    return result;
}

}  // namespace OpenABF
