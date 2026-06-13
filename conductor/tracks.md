# Tracks Registry

## Active Tracks

| Status | Track ID | Title | GitHub | Created | Updated |
| ------ | -------- | ----- | ------ | ------- | ------- |
| open | A7 | Multi-pin UV constraints for LSCM | [#42](https://github.com/educelab/OpenABF/issues/42) | 2026-03-16 | 2026-03-16 |
| open | P4 | IncompleteCholesky preconditioner for AngleBasedLSCM | [#47](https://github.com/educelab/OpenABF/issues/47) | 2026-03-18 | 2026-03-18 |
| open | P5 | HLSCM hot-path allocation reduction (UV map, vertexNeighbors, originalToLocal) | [#63](https://github.com/educelab/OpenABF/issues/63) | 2026-03-20 | 2026-03-20 |
| open | A8 | Extract shared LSCM system-building logic from solveLSCMLevel and ComputeImpl | [#64](https://github.com/educelab/OpenABF/issues/64) | 2026-03-20 | 2026-03-20 |
| open | P6 | Precompute sin/cos values before HLSCM face assembly loop | [#66](https://github.com/educelab/OpenABF/issues/66) | 2026-03-20 | 2026-03-20 |
| open | A10 | Add static Compute() overloads accepting levelRatio and minCoarseVertices | [#68](https://github.com/educelab/OpenABF/issues/68) | 2026-03-20 | 2026-03-20 |
| open | E2 | Add sphere cap built-in mesh generator to benchmark | [#48](https://github.com/educelab/OpenABF/issues/48) | 2026-03-20 | 2026-03-20 |
| open | F5 | Implement Hierarchical SLIM (H-SLIM) | [#52](https://github.com/educelab/OpenABF/issues/52) | 2026-03-20 | 2026-03-20 |
| open | F6 | Migrate to C++20 | [#54](https://github.com/educelab/OpenABF/issues/54) | 2026-03-20 | 2026-03-20 |
| open | F7 | Implement Policy-Based Math Backends | [#55](https://github.com/educelab/OpenABF/issues/55) | 2026-03-20 | 2026-03-20 |
| open | F8 | Implement CUDA/cuSolver Backend | [#56](https://github.com/educelab/OpenABF/issues/56) | 2026-03-20 | 2026-03-20 |
| open | F9 | Implement OpenBLAS/LAPACK Backend | [#57](https://github.com/educelab/OpenABF/issues/57) | 2026-03-20 | 2026-03-20 |
| open | F10 | Implement Ceres Solver Backend | [#58](https://github.com/educelab/OpenABF/issues/58) | 2026-03-20 | 2026-03-20 |
| open | F11 | Implement ACVD | [#62](https://github.com/educelab/OpenABF/issues/62) | 2026-03-20 | 2026-03-20 |
| open | F2 | Multi-chart UV packing | [#18](https://github.com/educelab/OpenABF/issues/18) | 2026-03-13 | 2026-03-13 |

## Archived Tracks

| Status | Track ID | Title | GitHub | Archived |
| ------ | -------- | ----- | ------ | -------- |
| closed | T5 | HLSCM internal component unit tests (buildHierarchy, prolongateUVs, solveLSCMLevel) | [#65](https://github.com/educelab/OpenABF/issues/65) | 2026-06-12 |
| closed | T7 | ABFPlusPlusAnglePreservation now uses curved mesh (hemisphere) for meaningful ABF exercise | [#50](https://github.com/educelab/OpenABF/issues/50) | 2026-06-12 |
| closed | A9 | Mark detail::hlscm types with @internal Doxygen tag | [#67](https://github.com/educelab/OpenABF/issues/67) | 2026-06-12 |
| closed | T6 | Add test for HLSCM::setPinnedVertices() instance API | [#51](https://github.com/educelab/OpenABF/issues/51) | 2026-06-12 |
| closed | T8 | Strengthen InstanceAPILevelRatio to verify parameters are applied | [#49](https://github.com/educelab/OpenABF/issues/49) | 2026-06-12 |
| closed | F3 | Multi-component extraction and parameterization pipeline (PRs #80, #81; ParameterizeConnectedComponents descoped) | [#19](https://github.com/educelab/OpenABF/issues/19) | 2026-06-12 |
| closed | M5 | split_edge and detail::filter cleanup (PR #79) | [#76](https://github.com/educelab/OpenABF/issues/76), [#77](https://github.com/educelab/OpenABF/issues/77), [#78](https://github.com/educelab/OpenABF/issues/78) | 2026-06-12 |
| closed | B1 | ABF lambda update index wrong for 2+ interior vertices | [#6](https://github.com/educelab/OpenABF/issues/6) | 2026-03-20 |
| closed | B2 | ABFPlusPlus LambdaStarInv inversion uses hardcoded `1.F` | [#7](https://github.com/educelab/OpenABF/issues/7) | 2026-03-20 |
| closed | B3 | Vec reverse iterators have wrong return type | [#8](https://github.com/educelab/OpenABF/issues/8) | 2026-03-20 |
| closed | B4 | Vec binary `*`/`/` inconsistent with `*=`/`/=` | [#9](https://github.com/educelab/OpenABF/issues/9) | 2026-03-20 |
| closed | B5 | `is_file_type` accesses `ext[0]` without bounds check | [#10](https://github.com/educelab/OpenABF/issues/10) | 2026-03-20 |
| closed | B6 | PLY reader vmap silently uninitialized if x/y/z absent | [#11](https://github.com/educelab/OpenABF/issues/11) | 2026-03-20 |
| closed | B7 | PlanGrad allocates and discards wheel vector | [#12](https://github.com/educelab/OpenABF/issues/12) | 2026-03-20 |
| closed | P1 | Mesh iteration methods return `std::vector` by value | [#4](https://github.com/educelab/OpenABF/issues/4) | 2026-03-20 |
| closed | P2 | `std::map` for interior-vertex index lookup in solver | [#13](https://github.com/educelab/OpenABF/issues/13) | 2026-03-20 |
| closed | P3 | `InitializeAnglesAndWeights` calls wheel() redundantly | [#14](https://github.com/educelab/OpenABF/issues/14) | 2026-03-20 |
| closed | A1 | Convergence tolerance hardcoded in ABF and ABFPlusPlus | [#15](https://github.com/educelab/OpenABF/issues/15) | 2026-03-20 |
| closed | A2 | FacePtr range iteration requires `*face` dereference | [#4](https://github.com/educelab/OpenABF/issues/4) | 2026-03-20 |
| closed | A3 | `insert_face` boundary update behavior undocumented | [#16](https://github.com/educelab/OpenABF/issues/16) | 2026-03-20 |
| closed | A4 | `gradient()` missing `[[nodiscard]]` | [#17](https://github.com/educelab/OpenABF/issues/17) | 2026-03-20 |
| closed | A5 | No configurable pinned edge selection for LSCM | [#2](https://github.com/educelab/OpenABF/issues/2) | 2026-03-20 |
| closed | A6 | AngleBasedLSCM A/B matrix construction hard to follow | [#3](https://github.com/educelab/OpenABF/issues/3) | 2026-03-20 |
| closed | F1 | HLSCM (Hierarchical LSCM) | [#5](https://github.com/educelab/OpenABF/issues/5) | 2026-03-20 |
| closed | F4 | Double-precision test coverage | [#20](https://github.com/educelab/OpenABF/issues/20) | 2026-03-20 |
| closed | T1 | ABF/ABFPlusPlus tests need multi-interior-vertex mesh | [#21](https://github.com/educelab/OpenABF/issues/21) | 2026-03-20 |
| closed | T2 | No test for FindEdgePath on disconnected mesh | [#22](https://github.com/educelab/OpenABF/issues/22) | 2026-03-20 |
| closed | T3 | No round-trip IO test for OBJ and PLY | [#23](https://github.com/educelab/OpenABF/issues/23) | 2026-03-20 |
| closed | T4 | No test for Vec reverse iterators | [#24](https://github.com/educelab/OpenABF/issues/24) | 2026-03-20 |
| closed | M1 | `Face::barycenter()` hardcodes `Vec<T, 3>` | [#25](https://github.com/educelab/OpenABF/issues/25) | 2026-03-20 |
| closed | M2 | `Edge::magnitude()` is non-const | [#26](https://github.com/educelab/OpenABF/issues/26) | 2026-03-20 |
| closed | M3 | `detail::erase_if` is dead code | [#27](https://github.com/educelab/OpenABF/issues/27) | 2026-03-20 |
| closed | M4 | `operator<<` for Vec defined outside OpenABF namespace | [#28](https://github.com/educelab/OpenABF/issues/28) | 2026-03-20 |
| closed | E1 | Flattening benchmark example | [#45](https://github.com/educelab/OpenABF/issues/45) | 2026-03-20 |
