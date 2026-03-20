# T5 Implementation Plan

## Phase 1: Expose level count
1. Change `buildHierarchy` return type to include the level count, or add a thin `buildHierarchyWithCount` wrapper
2. Update `MultiLevelHierarchy` test to assert `levels.size() >= 2` directly

## Phase 2: BuildHierarchy_LevelCount test
1. Call `detail::hlscm::buildHierarchy` on a 20×20 wavy surface with `minCoarseVerts=10`
2. Assert `levels.size() >= 3` (known from existing benchmark)
3. For each consecutive pair of levels, assert vertex count ratio ≈ `levelRatio`
4. For each level, verify `localToOriginal` and `originalToLocal` are consistent inverses
5. Assert pin vertices (idx 0 and 1) appear in every level's `localToOriginal`

## Phase 3: ProlongateUVs_BarycentricReconstruction test
1. Build a 2-level hierarchy on a small grid
2. Assign known UV positions to coarse-level vertices
3. Call `prolongateUVs` with the known UVs and the collapse records
4. For each collapsed vertex, compute expected UV from its `containingTri` + `bary` manually
5. Assert each prolongated UV matches expected within 1e-5

## Phase 4: SolveLSCMLevel_KnownMesh test
1. Build a single-level `HierarchyLevel` for the pyramid mesh
2. Call `solveLSCMLevel` with pins 0 and 1, no initial guess
3. Assert returned UV map has finite values and matches `AngleBasedLSCM::Compute` on same mesh
4. Assert pinned vertices have exactly the expected UV positions (same as LSCM pin selection)

## Phase 5: Verify and finalize
- `ctest` all tests pass
- `git clang-format`
