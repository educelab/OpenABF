#include <gtest/gtest.h>

#include "OpenABF/OpenABF.hpp"

using namespace OpenABF;

// The Parameterize helper needs ABF-traited edges/faces for the optimizer.
using Optimizer = ABFPlusPlus<float>;
using MeshType = Optimizer::Mesh;
using Param = AngleBasedLSCM<float, MeshType>;

TEST(HalfEdgeMeshUtils, ParameterizeWritesUVBackToOriginal)
{
    // Two disjoint triangles. ParameterizeConnectedComponents should flatten
    // each component and write UV coords back into the source's vertex
    // positions via the back-map.
    auto mesh = MeshType::New();
    mesh->insert_vertices({
        {0.f, 0.f, 0.f},
        {1.f, 0.f, 0.f},
        {0.f, 1.f, 0.f},
        {5.f, 5.f, 0.f},
        {6.f, 5.f, 0.f},
        {5.f, 6.f, 0.f},
    });
    mesh->insert_faces({{0, 1, 2}, {3, 4, 5}});

    std::vector<Vec<float, 3>> before;
    for (const auto& v : mesh->vertices()) {
        before.push_back(v->pos);
    }

    EXPECT_NO_THROW((ParameterizeConnectedComponents<Optimizer, Param>(mesh)));

    // Both components should have had at least one vertex rewritten.
    bool cc0Changed = false;
    bool cc1Changed = false;
    for (std::size_t i = 0; i < mesh->num_vertices(); ++i) {
        if (mesh->vertex(i)->pos != before[i]) {
            (i < 3 ? cc0Changed : cc1Changed) = true;
        }
    }
    EXPECT_TRUE(cc0Changed);
    EXPECT_TRUE(cc1Changed);
}
