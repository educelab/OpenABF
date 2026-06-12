/**
 * @example MultiChartFlatten.cpp
 *
 * # Multi-chart flattening demo
 *
 * Builds a mesh with multiple connected components (here, a 3x3 grid torn
 * along two seams), extracts each component as an independent mesh, runs
 * ABF++ + LSCM on each, and writes each flattened chart to its own .obj
 * file.
 *
 * The original mesh is never modified — the extracted sub-meshes own their
 * own vertices and per-edge state, and each is parameterized in isolation.
 *
 * @see OpenABF::HalfEdgeMesh::split_path
 * @see OpenABF::HalfEdgeMesh::extract_connected_components
 */
#include <iostream>
#include <string>

#include "OpenABF/OpenABF.hpp"

int main()
{
    using ABF = OpenABF::ABFPlusPlus<float>;
    using LSCM = OpenABF::AngleBasedLSCM<float, ABF::Mesh>;

    // Build a 3x3 grid (9 vertices, 8 triangles)
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

    // Tear the mesh along the vertical center line (v1 -> v4 -> v7) so the
    // grid splits into a left and right chart.
    std::cout << "Before split: " << mesh->num_connected_components() << " component(s)\n";
    mesh->split_path({1, 4});
    mesh->split_path({4, 7});
    std::cout << "After split:  " << mesh->num_connected_components() << " component(s)\n";

    // Extract each connected component as an independent mesh. The back-map
    // we get alongside each chart is the bridge between chart-vertex-idx and
    // the original-mesh vertex-idx; downstream code that needs to scatter
    // results back (e.g. a wedge UV table for mesh IO) uses it.
    auto charts = mesh->extract_connected_components();
    std::cout << "Extracted " << charts.size() << " chart(s)\n";

    // Flatten each chart in isolation and write it out as its own .obj.
    for (std::size_t i = 0; i < charts.size(); ++i) {
        auto& [chart, backMap] = charts[i];

        std::size_t iters{0};
        float grad{OpenABF::INF<float>};
        ABF::Compute(chart, iters, grad);
        LSCM::Compute(chart);

        const auto out = "openabf_example_multi_chart_" + std::to_string(i) + ".obj";
        OpenABF::WriteMesh(out, chart);
        std::cout << "Chart " << i << ": " << chart->num_vertices() << " vertices, "
                  << chart->num_faces() << " faces, " << iters << " ABF++ iters -> " << out << "\n";

        // back_map[chart_idx] -> original_idx, available for downstream uses
        // such as building a per-wedge UV map keyed by the source mesh.
        (void)backMap;
    }

    // The source mesh's 3D vertex positions are unchanged by the per-chart
    // flattening — only the extracted sub-meshes hold the 2D UV result.
    std::cout << "Source mesh 3D positions intact: " << mesh->num_vertices() << " vertices, "
              << mesh->num_faces() << " faces\n";
}
