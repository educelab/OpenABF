#pragma once

#include "OpenABF/HalfEdgeMesh.hpp"

namespace OpenABF
{

template <typename T, size_t Dims = 3>
class ABF
{
    struct VertProperties {
    };

    struct EdgeProperties {
        T alpha{0};
        T beta{0};
        T weight{0};
    };

    struct FaceProperties {
    };

    using HEM =
        HalfEdgeMesh<T, Dims, VertProperties, EdgeProperties, FaceProperties>;

public:
    HEM mesh;
};

}  // namespace OpenABF