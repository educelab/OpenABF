#pragma once

#include "OpenABF/OpenABF.hpp"

namespace OpenABF::tests
{
/** Construct a triangular pyramid with an open bottom */
template <typename MeshType>
auto ConstructPyramid() -> typename MeshType::Pointer
{
    auto mesh = MeshType::New();
    mesh->insert_vertex(0, 0, 0);
    mesh->insert_vertex(2, 0, 0);
    mesh->insert_vertex(1, std::sqrt(3), 0);
    mesh->insert_vertex(1, std::sqrt(3) / 3, std::sqrt(6) * 2 / 3);

    mesh->insert_faces({{1, 3, 0}, {3, 2, 0}, {3, 1, 2}});
    return mesh;
}
}  // namespace OpenABF::tests
