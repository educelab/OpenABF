# T5 Implementation Plan

## Phase 1: Expose level count — DONE (no header change required)
- `detail::hlscm::buildHierarchy` already returns the levels vector by value, so
  `levels.size()` directly yields the level count. No public-API change made.
- `MultiLevelHierarchy` test now calls `buildHierarchy` directly and asserts
  `levels.size() >= 2`.

## Phase 2: BuildHierarchy_LevelCount test — DONE
- Calls `detail::hlscm::buildHierarchy` on a 20×20 wavy surface with
  `levelRatio=4`, `minCoarseVerts=10`.
- Asserts `levels.size() >= 3`.
- Asserts vertex-count ratio between consecutive non-floor levels lies in
  `[0.5 * levelRatio, 2 * levelRatio]`.
- Asserts `localToOriginal` and `originalToLocal` are consistent inverses on
  every level.
- Asserts pins 0 and 1 appear in every level's `originalToLocal`.

## Phase 3: ProlongateUVs_BarycentricReconstruction test — DONE
- Builds a 2-level hierarchy on a 5×5 grid with `levelRatio=4`, `minCoarseVerts=5`.
- Assigns deterministic UVs (`U = origIdx, V = -origIdx`) to coarse vertices.
- Replays barycentric expansion in reverse-collapse order and asserts each
  prolongated UV matches within `1e-5f`.
- Also asserts coarse UVs survive untouched.

## Phase 4: SolveLSCMLevel_KnownMesh test — DONE
- Builds a single-level pyramid `HierarchyLevel`.
- Calls `detail::hlscm::solveLSCMLevel<float, ConjugateGradient>` with pins 0/1
  and no initial guess.
- Asserts finite UVs, pin0 at `(0, 0)`, pin1 at `(2, 0)`.
- Asserts result matches `AngleBasedLSCM<float, ..., ConjugateGradient>::Compute`
  within `1e-4f` on the same mesh.

## Phase 5: Verify and finalize — DONE
- `ctest`: all 6 test executables pass.
- `git clang-format` applied to `tests/src/TestParameterization.cpp`.
- No `include/` changes → no single-header regeneration required.
