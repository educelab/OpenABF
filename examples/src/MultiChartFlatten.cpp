/**
 * @example MultiChartFlatten.cpp
 *
 * # Multi-chart flattening demo
 *
 * Builds a mesh with multiple connected components (here, a 3x3 grid torn
 * down its center line), extracts each component as an independent mesh, runs
 * ABF++ + LSCM on each, packs the flattened charts into a shared coordinate
 * frame with OpenABF::PackCharts, and writes the packed atlas to a single
 * .obj file.
 *
 * Charts are flattened defensively: not every chart topology is solvable, so a
 * chart that throws OpenABF::SolverException is reported and skipped rather
 * than aborting the atlas.
 *
 * The original mesh is never modified — the extracted sub-meshes own their
 * own vertices and per-edge state, and each is parameterized in isolation
 * before being placed into the common frame.
 *
 * @see OpenABF::HalfEdgeMesh::split_path
 * @see OpenABF::HalfEdgeMesh::extract_connected_components
 * @see OpenABF::PackCharts
 */
#include <iostream>
#include <vector>

#include "OpenABF/OpenABF.hpp"

int main()
{
    using ABF = OpenABF::ABFPlusPlus<float>;
    using LSCM = OpenABF::AngleBasedLSCM<float, ABF::Mesh>;
    using Mesh = ABF::Mesh;

    // Build a 3x3 grid (9 vertices, 8 triangles)
    auto mesh = Mesh::New();
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

    // Tear the mesh along the vertical center line (v1 -> v4 -> v7) so the
    // grid splits into a left and right chart.
    std::cout << "Before split: " << mesh->num_connected_components() << " component(s)\n";
    mesh->split_path({1, 4, 7});
    std::cout << "After split:  " << mesh->num_connected_components() << " component(s)\n";

    // Extract each connected component as an independent mesh. Each chart
    // comes with a `vertex_map` and `face_map` that bridge chart indices to
    // source-mesh indices; downstream code that scatters per-vertex or
    // per-face data back to the source (e.g. wedge UV tables, per-face
    // material assignments) uses them.
    auto charts = mesh->extract_connected_components();
    std::cout << "Extracted " << charts.size() << " chart(s)\n";

    // Flatten each chart in isolation, collecting the parameterized meshes.
    std::vector<Mesh::Pointer> chartMeshes;
    for (std::size_t i = 0; i < charts.size(); ++i) {
        auto& cc = charts[i];

        std::size_t iters{0};
        float grad{OpenABF::INF<float>};
        // Not every chart topology is solvable, and one unsolvable chart should
        // not cost you the whole atlas. Flatten each chart inside its own
        // try/catch and pack whatever succeeded: a chart that throws is simply
        // left out of `chartMeshes`, so its faces get no UVs while the rest of
        // the atlas is still produced.
        try {
            ABF::Compute(cc.mesh, iters, grad);
            LSCM::Compute(cc.mesh);
            std::cout << "Chart " << i << ": " << cc.mesh->num_vertices() << " vertices, "
                      << cc.mesh->num_faces() << " faces, " << iters << " ABF++ iters\n";
            chartMeshes.push_back(cc.mesh);
        } catch (const OpenABF::SolverException& e) {
            std::cout << "Chart " << i << ": skipped, could not be flattened (" << e.what()
                      << ")\n";
        }

        // cc.vertex_map[chart_idx] -> original vertex idx
        // cc.face_map[chart_idx]   -> original face idx
        // Available for downstream uses such as building a per-wedge UV map
        // keyed by source-mesh face corners.
    }

    // Pack the flattened charts into a shared frame. `normalize` fits the whole
    // atlas into [0,1]^2 via a single global uniform scale, which preserves the
    // charts' relative sizes. Packing only edits each chart's 2D vertex
    // positions in place; the per-chart vertex_map/face_map remain valid.
    OpenABF::PackOptions<float> opts;
    opts.normalize = true;
    // Add a gutter around every chart so neighbouring charts -- and the atlas
    // boundary -- don't touch. Without this, charts pack flush and texture
    // filtering can bleed one chart's texels into another. `padding` is in
    // absolute chart units and is applied before the normalize scaling.
    opts.padding = 0.1f;
    auto extent = OpenABF::PackCharts<Mesh>(chartMeshes, opts);
    std::cout << "Packed atlas extent: [" << extent.min[0] << ", " << extent.min[1] << "] -> ["
              << extent.max[0] << ", " << extent.max[1] << "]\n";

    // Merge the packed charts into a single mesh and write it as one atlas.
    // MergeMeshes returns provenance maps (vertex_source/face_source) that, when
    // composed with each component's vertex_map/face_map, trace any atlas
    // element back to the torn source mesh.
    auto merged = OpenABF::MergeMeshes<Mesh>(chartMeshes);

    const std::string out = "openabf_example_multi_chart_packed.obj";
    OpenABF::WriteMesh(out, merged.mesh);
    std::cout << "Wrote packed atlas: " << merged.mesh->num_vertices() << " vertices, "
              << merged.mesh->num_faces() << " faces -> " << out << "\n";

    /*
     * Reference: building a per-wedge UVMap from the merged result
     * -----------------------------------------------------------------------
     * OpenABF does not own a UV-map type, but the merged atlas plus the
     * components' back-maps carry everything needed to populate one. The
     * snippet below (not compiled here) targets educelab::core's UVMap:
     *
     *     educelab/core/types/UVMap.hpp
     *
     * The UVMap is keyed by (face, corner) against the *torn source mesh*
     * `mesh` — which still holds the original 3D geometry, since only the
     * extracted charts were flattened. UV coordinates come from the packed
     * chart vertices. Corner positions are resolved by *vertex identity*, not
     * by traversal order: a face's winding may be reversed at insertion time,
     * so the chart/atlas corner order is not guaranteed to match the source
     * face's corner order (see PackCharts / HalfEdgeMesh::insert_face).
     *
     *     #include <algorithm>
     *     #include "educelab/core/types/UVMap.hpp"
     *     using educelab::UVMap;
     *
     *     UVMap<float, 2> uv;
     *
     *     for (std::size_t mf = 0; mf < merged.mesh->num_faces(); ++mf) {
     *         // Atlas face -> source chart + chart-local face -> source (M') face.
     *         const auto [chart, subFace] = merged.face_source[mf];
     *         const auto srcFace = charts[chart].face_map[subFace];
     *
     *         // Source face corner order, keyed by M' vertex index.
     *         std::vector<std::size_t> srcCorners;
     *         for (const auto& e : *mesh->faces()[srcFace]) {
     *             srcCorners.push_back(e->vertex->idx);
     *         }
     *
     *         // Each atlas-face corner carries its packed UV in pos.
     *         for (const auto& e : *merged.mesh->faces()[mf]) {
     *             const auto [vChart, vSub] = merged.vertex_source[e->vertex->idx];
     *             const auto srcVert = charts[vChart].vertex_map[vSub];  // M' vertex
     *
     *             // Place the UV at the matching corner of the source face.
     *             const auto corner = static_cast<std::size_t>(std::distance(
     *                 srcCorners.begin(),
     *                 std::find(srcCorners.begin(), srcCorners.end(), srcVert)));
     *
     *             uv.map(srcFace, corner, uv.insert(e->vertex->pos[0],
     *                                               e->vertex->pos[1]));
     *         }
     *     }
     *
     *     // uv.get_coordinate(srcFace, corner) now yields the packed UV for
     *     // each wedge of `mesh`, ready for OBJ `vt` emission.
     *
     * This table is valid for BOTH the torn mesh and the untorn (pre-split)
     * mesh. split_path preserves face indices and per-face winding (it never
     * re-inserts faces), and the (face, corner) keys are resolved by vertex
     * identity against `mesh`'s own faces — so any winding reversal insert_face
     * applies (when M is built, when components are cloned by extract, and when
     * charts are merged) is absorbed rather than baked into the keys. The rule
     * that makes this work: consume the UVMap against the same HalfEdgeMesh (or
     * one sharing its winding) and resolve corners by identity, never by a raw
     * traversal index. Caveat: insert_face's auto-rewinding means a face's
     * as-built corner order may differ from the raw input face list, and the
     * mesh does not record that permutation (see issue tracker / bug track B9).
     */

    // The source mesh's 3D vertex positions are unchanged by the per-chart
    // flattening and packing — only the extracted sub-meshes hold the 2D UV
    // result.
    std::cout << "Source mesh 3D positions intact: " << mesh->num_vertices() << " vertices, "
              << mesh->num_faces() << " faces\n";
}
