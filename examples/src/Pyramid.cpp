#include <iostream>

#include "OpenABF/OpenABF.hpp"

int main()
{
    using ABF = OpenABF::ABFPlusPlus<float>;
    using LSCM = OpenABF::AngleBasedLSCM<float, ABF::Mesh>;

    auto mesh = ABF::Mesh::New();
    mesh->insert_vertex(0, 0, 0);
    mesh->insert_vertex(2, 0, 0);
    mesh->insert_vertex(1, std::sqrt(3), 0);
    mesh->insert_vertex(1, std::sqrt(3) / 3, std::sqrt(6) * 2 / 3);

    mesh->insert_face(0, 3, 1);
    mesh->insert_face(0, 2, 3);
    mesh->insert_face(2, 1, 3);

    std::vector<OpenABF::Vec3f> orig;
    for (const auto& v : mesh->vertices()) {
        orig.push_back(v->pos);
    }

    std::size_t iters{0};
    float grad{OpenABF::INF<float>};
    ABF::Compute(mesh, iters, grad);

    std::cout << "ABF++ Final gradient: " << grad << std::endl;
    std::cout << "ABF++ Iterations: " << iters << std::endl;

    LSCM::Compute(mesh);

    std::vector<OpenABF::Vec3f> uvs;
    for (const auto& v : mesh->vertices()) {
        std::cout << v->idx << ": " << orig[v->idx] << " -> " << v->pos
                  << std::endl;
    }
}
