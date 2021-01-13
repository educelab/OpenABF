#include <iostream>

#include "OpenABF/ABF.hpp"
#include "OpenABF/HalfEdgeMesh.hpp"

int main()
{
    OpenABF::ABF<float> abf;

    abf.mesh.insertVertex(0, 0, 0);
    abf.mesh.insertVertex(1, 0, 0);
    abf.mesh.insertVertex(1, -1, 0);
    abf.mesh.insertVertex(1, 1, 0);

    abf.mesh.insertFace(0, 1, 2);
    abf.mesh.insertFace(0, 3, 1);
    abf.mesh.insertFace(0, 2, 3);

    std::cout << "Interior vertices: [edges]" << std::endl;
    size_t count{0};
    for (const auto& v : abf.mesh.interior_vertices()) {
        std::cout << "\t- Vertex " << v->idx << ": ";
        for (const auto& e : v->edges()) {
            if (count++ > 0) {
                std::cout << ", ";
            }
            std::cout << e->vertex->idx << "<->" << e->next->vertex->idx;
        }
    }
    std::cout << std::endl;

    std::cout << "Boundary vertices: ";
    count = 0;
    for (const auto& v : abf.mesh.boundary_vertices()) {
        if (count++ > 0) {
            std::cout << ", ";
        }
        std::cout << v->idx;
    }
    std::cout << std::endl;
}
