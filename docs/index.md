**OpenABF** is a header-only C++ library of angle-based flattening algorithms.
It is designed to be as simple as possible to integrate into existing projects.

## Requirements
* [Eigen 3.3+](http://eigen.tuxfamily.org/)

## Installation
### CMake
This project is built and installed using the CMake build system:

```{.sh}
mkdir build
cmake -S . -B build/
cmake --install build/
```

This will install the OpenABF header(s) to your default system path and provide
an easy method to link OpenABF against your own CMake project:

```cmake
# Find OpenABF libraries
find_package(OpenABF REQUIRED)

# Link to an executable
add_executable(MyTarget main.cpp)
target_link_libraries(MyTarget OpenABF::OpenABF)
```

The CMake project provides a number of flags for configuring the installation:
- `OPENABF_MULTIHEADER`: Install the multi-header version of OpenABF
  (Default: OFF)
- `OPENABF_BUILD_TESTS`: Build project unit tests. This will download and build
  the Google Test framework. (Default: OFF)
- `OPENABF_BUILD_DOCS`: Build documentation. Dependencies: Doxygen, Graphviz
  (optional). (Default: ON if Doxygen is found)

### Manual
Copy and paste the contents of `single_include` to your project or include path.
As this project requires Eigen, you also need to add that project to your
include path.

## Usage
```{.cpp}
#include <OpenABF/OpenABF.hpp>

// Alias algorithms for convenience
using ABF = OpenABF::ABFPlusPlus<float>;
using LSCM = OpenABF::AngleBasedLSCM<float, ABF::Mesh>;

// Make a new mesh
auto mesh = ABF::Mesh::New();
mesh->insert_vertex(0, 0, 0);
mesh->insert_vertex(2, 0, 0);
mesh->insert_vertex(1, std::sqrt(3), 0);
mesh->insert_vertex(1, std::sqrt(3) / 3, 1);

mesh->insert_face(0, 3, 1);
mesh->insert_face(0, 2, 3);
mesh->insert_face(2, 1, 3);

// Print original coordinates
for (const auto& v : mesh->vertices()) {
    std::cout << v->idx << ": " << v->pos << std::endl;
}

// Compute parameterized angles
ABF abf;
abf.setMesh(mesh);
abf.compute();

// Compute mesh parameterization from angles
LSCM lscm;
lscm.setMesh(mesh);
lscm.compute();

// Print new coordinates
for (const auto& v : mesh->vertices()) {
    std::cout << v->idx << ": " << v->pos << std::endl;
}
```