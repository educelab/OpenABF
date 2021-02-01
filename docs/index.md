**%OpenABF** is a single-header C++ library of angle-based flattening algorithms.
The templated interface is designed for simple out-of-the-box use, and 
integration with existing geometric processing pipelines is quick and easy.
View the latest source code on [GitLab](https://gitlab.com/educelab/OpenABF).

## Dependencies
- C++14 compiler
- [Eigen 3.3+](http://eigen.tuxfamily.org/)
- CMake 3.15+ (optional)

## Usage
The following example demonstrates how to construct and parameterize a mesh 
with OpenABF: 
 
```{.cpp}
#include <OpenABF/OpenABF.hpp>

// Alias algorithms for convenience
using ABF = OpenABF::ABFPlusPlus<float>;
using LSCM = OpenABF::AngleBasedLSCM<float, ABF::Mesh>;

// Make a triangular pyramid mesh
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
ABF::Compute(mesh);

// Compute mesh parameterization from angles
LSCM::Compute(mesh);

// Print new coordinates
for (const auto& v : mesh->vertices()) {
    std::cout << v->idx << ": " << v->pos << std::endl;
}
```
**Note:** The `HalfEdgeMesh` class 
[currently assumes](https://gitlab.com/educelab/OpenABF/-/issues/4) that the 
surface has a boundary, is manifold, and that the winding order of all faces is 
the same. Care should be taken that this assumption is not violated when 
constructing your mesh.

## Installation
### CMake
This project can be configured and installed using the CMake build system:

```{.sh}
mkdir build
cmake -S . -B build/
cmake --install build/
```

This will install the OpenABF header(s) to your system include path and provide 
an easy method for including OpenABF inside of your own CMake project:

```
# Find OpenABF libraries
find_package(OpenABF REQUIRED)

# Link to an executable
add_executable(MyTarget main.cpp)
target_link_libraries(MyTarget OpenABF::OpenABF)
```
  
**Note:** For best performance, configure your CMake project with the 
`-DCMAKE_BUILD_TYPE=Release` flag.

#### Configuration
The OpenABF CMake project provides a number of flags for configuring the 
installation:
- `OPENABF_MULTIHEADER`: Install the multi-header version of OpenABF 
  (Default: OFF)
- `OPENABF_BUILD_EXAMPLES`: Build example applications. (Default: OFF)
- `OPENABF_BUILD_TESTS`: Build project unit tests. This will download and build
  the Google Test framework. (Default: OFF)
- `OPENABF_BUILD_DOCS`: Build documentation. Dependencies: Doxygen, Graphviz
  (optional). Unavailable if Doxygen is not found. (Default: OFF)

#### FetchContent (CMake 3.11+)
Another option for providing OpenABF to your project is by using CMake's 
[FetchContent module](https://cmake.org/cmake/help/latest/module/FetchContent.html):

```
include(FetchContent)
FetchContent_Declare(
  openabf
  GIT_REPOSITORY https://gitlab.com/educelab/OpenABF.git
  GIT_TAG v1.0
)

# Populate the project but exclude from All targets
FetchContent_GetProperties(openabf)
if(NOT openabf_POPULATED)
  FetchContent_Populate(openabf)
  add_subdirectory(${openabf_SOURCE_DIR} ${openabf_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()
```

This downloads the OpenABF source code and adds it to your CMake project as a 
subproject. Link it against your targets as you would any library added with 
`find_package`:

```
add_executable(MyTarget main.cpp)
target_link_libraries(MyTarget OpenABF::OpenABF)
```

### Manual
Copy and paste the contents of `single_include/` to your project or include 
path. As OpenABF depends upon the Eigen library, you will also need to add the 
Eigen headers to your include path:

```{.sh}
g++ -I /path/to/eigen/ -std=c++14 -DNDEBUG -O3 main.cpp -o main
```

**Note:** For best performance, compile your application with the `-DNDEBUG -03` 
preprocessor definitions.