#include <iostream>

#include "abf++/ABF.hpp"
#include "abf++/HalfEdgeMesh.hpp"

int main()
{
    HalfEdgeMesh<> hem;

    hem.insertVertex(0, 0, 0);
    hem.insertVertex(1, 0, 0);
    hem.insertVertex(1, -1, 0);
    hem.insertVertex(1, 1, 0);

    hem.insertFace(0, 1, 2);
    hem.insertFace(0, 3, 1);

    std::cout << "End." << std::endl;
}
