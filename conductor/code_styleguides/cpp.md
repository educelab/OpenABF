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
- **Types / classes**: `PascalCase`
- **Functions / methods**: `camelCase` or `snake_case` — be consistent within a class
- **Member variables**: `camelCase` or `snake_case_` with trailing underscore for private members
- **Constants / enums**: `kPascalCase` or `ALL_CAPS`
- **Namespaces**: `lowercase`
- **Template parameters**: `PascalCase` (e.g. `typename Scalar`)

## C++ Standards
- **Standard**: C++17 required; do not use C++20 features in public headers without a guard
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
- **Formatter**: `clang-format` (use project root `.clang-format`)
- **Static analysis**: `clang-tidy` recommended (no project `.clang-tidy` yet — add one as needed)
- **CI**: GitHub Actions runs lint checks; all PRs must pass before merge
