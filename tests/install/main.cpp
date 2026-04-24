#include <cstdlib>
#include <iostream>

#include <OpenABF/OpenABF.hpp>

int main()
{
    // Verify OpenABF types work (depends on Eigen being found transitively)
    OpenABF::Vec3f v{1.0f, 2.0f, 3.0f};
    if (v[0] != 1.0f || v[1] != 2.0f || v[2] != 3.0f) {
        std::cerr << "Vec3f value mismatch" << std::endl;
        return EXIT_FAILURE;
    }
    std::cout << "OpenABF find_package test passed" << std::endl;
    return EXIT_SUCCESS;
}
