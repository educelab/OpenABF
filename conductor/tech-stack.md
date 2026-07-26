# Tech Stack

## Primary Language
- **C++20** — all library code; required standard for consumers
- **Python 3.x** — utility/build scripts (e.g. `thirdparty/amalgamate/amalgamate.py` for single-header generation); Python bindings (pybind11) are a future consideration

## Build System
- **CMake 3.15+** — primary build system; supports FetchContent, GNUInstallDirs, and optional subdirectories (docs, examples, tests)

## Library Type
- **Header-only** — distributed as a single amalgamated header (`single_include/OpenABF/OpenABF.hpp`) or multi-header form (`include/OpenABF/`)

## Dependencies
- **Eigen3** (≥3.3, prefers ≥5) — linear algebra; required at configure time
- **Google Test** — testing only, fetched via CMake FetchContent
- **Doxygen** (optional) — documentation generation

## Frontend / Backend / Database
- None — this is a C++ library with no runtime services

## Distribution
- **CMake FetchContent / manual copy** — primary integration method for consumers
- **Package managers** — vcpkg, Conan, and/or system packages (e.g. apt)

## Optional Future Additions
- Python bindings via pybind11
