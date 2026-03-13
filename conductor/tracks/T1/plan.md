# T1 Implementation Plan

## Phase 1 — Add grid fixture
1. In `tests/src/Utils.hpp`, add `ConstructGrid<MeshType>(rows, cols)` that builds
   a flat rectangular triangulated mesh
2. A 3×3 vertex grid (2×2 quads, 8 triangles) gives 1 interior vertex; a 4×4 grid
   (3×3 quads, 18 triangles) gives 4 interior vertices — use the latter

## Phase 2 — Add tests
1. `TEST(Parameterizations, ABF_MultiInterior)` — run ABF + LSCM on the grid, verify
   UV coordinates are plausible (no NaN, boundary stays on border)
2. `TEST(Parameterizations, ABFPlusPlus_MultiInterior)` — same
3. `TEST(Parameterizations, ABF_NoInteriorVertices)` — single triangle, no interior
   vertices; verify no crash and boundary vertices are placed correctly
4. `TEST(Parameterizations, ABF_MaxIters)` — set maxIters=1, verify that the solver
   exits after 1 iteration and `iterations() == 1`

## Phase 3 — Verify
- All new tests pass (may initially fail before B1 is fixed — that is expected and desired)
- Run `ctest` to confirm no regressions
