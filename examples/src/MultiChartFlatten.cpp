/**
 * @example MultiChartFlatten.cpp
 *
 * # Multi-chart flattening demo
 *
 * Builds a mesh with multiple connected components (here, a 3x3 grid torn
 * along two seams), extracts each component as an independent mesh, runs
 * ABF++ + LSCM on each, packs the flattened charts into a shared coordinate
 * frame with OpenABF::PackCharts, and writes the packed atlas to a single
 * .obj file.
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
        ABF::Compute(cc.mesh, iters, grad);
        LSCM::Compute(cc.mesh);
        chartMeshes.push_back(cc.mesh);

        std::cout << "Chart " << i << ": " << cc.mesh->num_vertices() << " vertices, "
                  << cc.mesh->num_faces() << " faces, " << iters << " ABF++ iters\n";

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
    auto extent = OpenABF::PackCharts<Mesh>(chartMeshes, opts);
    std::cout << "Packed atlas extent: [" << extent.min[0] << ", " << extent.min[1] << "] -> ["
              << extent.max[0] << ", " << extent.max[1] << "]\n";

    // Merge the packed charts into a single mesh and write it as one atlas.
    // The charts are disjoint in the packed frame, so concatenating their
    // vertices and re-emitting their faces (with an index offset per chart)
    // reconstructs a multi-component mesh in UV space.
    auto packed = Mesh::New();
    std::size_t offset{0};
    for (const auto& chart : chartMeshes) {
        for (const auto& v : chart->vertices()) {
            packed->insert_vertex(v->pos[0], v->pos[1], v->pos[2]);
        }
        for (const auto& face : chart->faces()) {
            std::vector<std::size_t> idxs;
            for (const auto& edge : *face) {
                idxs.push_back(edge->vertex->idx + offset);
            }
            packed->insert_face(idxs);
        }
        offset += chart->num_vertices();
    }

    const std::string out = "openabf_example_multi_chart_packed.obj";
    OpenABF::WriteMesh(out, packed);
    std::cout << "Wrote packed atlas: " << packed->num_vertices() << " vertices, "
              << packed->num_faces() << " faces -> " << out << "\n";

    // The source mesh's 3D vertex positions are unchanged by the per-chart
    // flattening and packing — only the extracted sub-meshes hold the 2D UV
    // result.
    std::cout << "Source mesh 3D positions intact: " << mesh->num_vertices() << " vertices, "
              << mesh->num_faces() << " faces\n";
}
