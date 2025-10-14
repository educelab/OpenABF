/**
 * @example FlattenTool.cpp
 *
 * # Flattening tool
 *
 * Demo utility for flattening the geometry of an OBJ mesh.
 */
#include <filesystem>

#include "OpenABF/OpenABF.hpp"

using T = float;
using ABF = OpenABF::ABFPlusPlus<T>;
using LSCM = OpenABF::AngleBasedLSCM<T, ABF::Mesh>;

namespace fs = std::filesystem;

auto main(const int argc, char* argv[]) -> int
{
    if (argc < 3) {
        std::cerr << "Usage: " << fs::path(argv[0]).filename().c_str()
                  << " [input.(obj|ply)] [output.(obj|ply)]";
        std::cerr << "\n";
        return EXIT_FAILURE;
    }

    const fs::path input_path = argv[1];
    const fs::path output_path = argv[2];

    std::cout << "Reading mesh...\n";
    auto mesh = OpenABF::ReadMesh<ABF::Mesh>(input_path);

    std::cout << "Mesh info:\n";
    std::cout << " - Num. vertices: " << mesh->num_vertices() << "\n";
    std::cout << " - Num. faces: " << mesh->num_faces() << "\n\n";

    std::cout << "Computing ABF++...\n";
    std::size_t iters{0};
    T grad{OpenABF::INF<T>};
    ABF::Compute(mesh, iters, grad);
    std::cout << "ABF++ :: Iterations: " << iters << " :: Gradient: " << grad
              << "\n";

    std::cout << "Computing LSCM...\n";
    LSCM::Compute(mesh);

    std::cout << "Writing mesh...\n";
    OpenABF::WriteMesh(output_path, mesh);
    std::cout << "Done.\n";
}