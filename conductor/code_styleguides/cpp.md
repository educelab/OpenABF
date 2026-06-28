# C++ Style Guide

Enforced via `.clang-format` at the project root. Run `clang-format` before committing.
The `thirdparty/` directory is excluded from formatting — it contains external code.

## Formatting (from .clang-format)
- **Base style**: Google (modified)
- **Indent width**: 4 spaces (no tabs)
- **Column limit**: 100 characters
- **Brace style**: Linux (opening brace on same line for functions, new line for control flow)
- **Pointer alignment**: Left (`int* p`, not `int *p`)
- **Template declarations**: Always break before (`template <typename T>` on its own line)
- **Short functions**: Allowed on a single line
- **Short if/for/while**: Not allowed on a single line
- **Include sorting**: Enabled; priority order: standard library, C headers, third-party, project headers

## Naming Conventions

| Identifier kind | Style | Examples |
|---|---|---|
| Types, classes, structs, type aliases, enum classes | `PascalCase` | `HalfEdgeMesh`, `PinMap`, `SolverException`, `DefaultVertexTraits` |
| Static class methods (factories, public algorithms exposed on a class) | `PascalCase` | `HalfEdgeMesh::New`, `ABF::Compute`, `MeshIO::Read` |
| Instance methods | `lower_snake_case` | `insert_vertex`, `is_boundary`, `set_pins`, `compute` |
| Free-function **algorithms** (named computational steps) | `PascalCase` | `ComputeMeshAngles`, `FindEdgePath`, `IsManifold` |
| Free-function **utilities** (small helpers, string/container ops) | `lower_snake_case` | `vec_to_string`, `remove_if`, `trim`, `icase_compare` |
| `detail::` namespace types | `PascalCase` | `detail::lscm::PinMap`, `detail::hlscm::DecimationMesh` |
| `detail::` namespace free functions | Follow the algorithm/utility split above | `detail::lscm::BuildSystem` (algorithm), `detail::trim` (utility) |
| `detail::` internal class methods | `lower_snake_case` | `DecimationMesh::try_collapse`, `DecimationMesh::is_alive` |
| Constants and enum values | `kPascalCase`, `ALL_CAPS`, or lowercase template constants — consistent within an enum/namespace | `PI<T>`, `Norm::L1`, `kMaxIters` |
| Namespaces | `lowercase` | `OpenABF`, `detail`, `lscm`, `hlscm`, `traits` |
| Template parameters | `PascalCase` | `T`, `MeshType`, `Solver` |
| Public data members (POD / trait types) | `lower_snake_case` | `Vertex::pos`, `Edge::alpha`, `Face::head` |
| Private and protected members | `lower_snake_case_` (trailing underscore) | `verts_`, `pinned_indices_`, `is_pinned_` |

### Algorithm vs utility — heuristic
A free function is an **algorithm** if a user would plausibly call it as a meaningful step in a flattening pipeline (or other domain workflow). It's a **utility** if it's a small generic helper that supports algorithms but isn't itself one. When in doubt: PascalCase for things that take a mesh/face/edge or compute something domain-specific; snake_case for things that take a string or generic container.

### Setter style
Setters on instance configuration follow the instance-method rule: `set_<thing>(value)`, not `setThing(value)`. Example: `set_level_ratio(4)`, not `setLevelRatio(4)`.

### Project history note (2026-06)
Pre-2026 code in several places uses `camelCase` for instance methods and `detail::` free functions (e.g., `setLevelRatio`, `detail::lscm::buildSystem`, `DecimationMesh::tryCollapse`). These are tracked for rename in issue #94. New code must follow the table above.

## C++ Standards
- **Standard**: C++20 required
- **Headers**: All public headers must be self-contained (include what they use)
- **`#pragma once`**: Preferred over include guards in project headers
- **RAII**: Prefer RAII over manual resource management
- **`auto`**: Use where it improves readability; avoid when the type is non-obvious
- **Exceptions**: Permitted for error signaling in the public API; use project exception types from `Exceptions.hpp`

## Library-Specific Rules
- Do not introduce runtime dependencies beyond Eigen3 in public headers
- Algorithms must be templated on scalar type (e.g. `float`, `double`) where applicable
- Internal implementation details belong in `detail/` namespaces or anonymous namespaces
- Avoid `using namespace` in headers
- All public API must have Doxygen-compatible documentation comments (`/** ... */`)

## Tooling
- **Formatter**: `clang-format` (use project root `.clang-format`). **Note**: CI installs `clang-format` via `apt` on `ubuntu-latest`, which currently resolves to v18; locally-installed newer versions (v20+) may produce different reflow on edge cases. If CI fails formatting after a local pre-commit format, manually match CI by inspecting the diff in the failed action log.
- **Static analysis**: `clang-tidy` recommended (no project `.clang-tidy` yet — add one as needed)
- **CI**: GitHub Actions runs lint checks; all PRs must pass before merge
