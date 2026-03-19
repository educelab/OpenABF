# Installation        {#installation}

%OpenABF can be installed through CMake or by manually copying the header 
file(s) into your include directory. This page outlines the various installation 
options.

## Dependencies
- C++17 compiler
- [Eigen 3.3+](http://eigen.tuxfamily.org/)
- CMake 3.15+ (optional)

## CMake
This project can be configured and installed using the CMake build system:

```{.sh}
mkdir build
cmake -S . -B build/
cmake --install build/
```

This will install the %OpenABF header(s) to your system include path and provide
an easy method for including OpenABF inside your own CMake project:

```
# Find OpenABF libraries
find_package(OpenABF REQUIRED)

# Link to an executable
add_executable(MyTarget main.cpp)
target_link_libraries(MyTarget OpenABF::OpenABF)
```

**Note:** For best performance, configure your CMake project with the
`-DCMAKE_BUILD_TYPE=Release` flag.

### Configuration
The %OpenABF CMake project provides a number of flags for configuring the
installation:
- `OPENABF_MULTIHEADER`: Install the multi-header version of %OpenABF
  (Default: OFF)
- `OPENABF_BUILD_EXAMPLES`: Build example applications. (Default: OFF)
- `OPENABF_BUILD_TESTS`: Build project unit tests. This will download and build
  the Google Test framework. (Default: OFF)
- `OPENABF_BUILD_DOCS`: Build documentation. Dependencies: Doxygen, Graphviz
  (optional). Unavailable if Doxygen is not found. (Default: OFF)

### FetchContent (CMake 3.11+)
Another option for providing OpenABF to your project is by using CMake's
[FetchContent module](https://cmake.org/cmake/help/latest/module/FetchContent.html):

```
include(FetchContent)
FetchContent_Declare(
  openabf
  GIT_REPOSITORY https://gitlab.com/educelab/OpenABF.git
  GIT_TAG v2.1.0
  EXCLUDE_FROM_ALL
)
FetchContent_MakeAvailable()
```

This downloads the %OpenABF source code and adds it to your CMake project as a
subproject. Link it against your targets as you would any library added with
`find_package`:

```
add_executable(MyTarget main.cpp)
target_link_libraries(MyTarget OpenABF::OpenABF)
```

## Manual installation (i.e. copy-and-paste)
Copy and paste the contents of `single_include/` to your project or include
path. As %OpenABF depends upon the Eigen library, you will also need to add the
Eigen headers to your include path:

```{.sh}
g++ -I /path/to/eigen/ -std=c++17 -DNDEBUG -O3 main.cpp -o main
```

**Note:** For best performance, compile your application with the `-DNDEBUG -03`
preprocessor definitions.

## Compilation on Windows

For many legacy reasons, the Microsoft Visual C++ compiler (MSVC) is not
automatically conformant with the C++ standard in all cases. This may lead to
the following issues when compiling against OpenABF.

**Note:** As this project only supports C++17 and up, you should always compile
with at least `/std:c++17`.

### Undeclared identifier errors

At the time of this writing, MSVC does not automatically recognize the
alternative boolean operators `and`, `or`, and `not`. There are multiple ways to
fix this issue.

- **Compile with the `/permissive-` flag:** This flag enables C++ language
  conformance for the entire compilation. This is currently Microsoft's suggested
  method for enabling the alternative operators, and is added by default when
  creating new projects in Visual Studio 2017 version 15.5 and later (but not when
  using MSVC on the command line). Since the `/permissive-` flag enables strict
  C++ conformance for the entire project (it enables much more than just operator
  support), it may lead to other compilation problems in existing code.
  [See the flag documentation](https://learn.microsoft.com/en-us/cpp/build/reference/permissive-standards-conformance)
  to decide if this is the right solution for your project.

- **Include `iso646.h`:** This header file provides definitions for the
  alternative operators and should be included before including OpenABF. This
  was Microsoft's previous recommendation for enabling alternative operator
  support. It may have fewer side effects in existing code bases than
  the `/permissive-` flag.
