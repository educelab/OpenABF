# T5 — HLSCM internal component unit tests

## GitHub Issue
https://github.com/educelab/OpenABF/issues/65

## Problem
All HLSCM tests currently go through the top-level `Compute()` API. Three internal
components have complex invariants that are difficult to diagnose via integration tests alone:

- `buildHierarchy` — level count, vertex counts, `localToOriginal`/`originalToLocal` consistency
- `prolongateUVs` — barycentric UV reconstruction per collapsed vertex
- `solveLSCMLevel` — LSCM solve on a small known mesh; pinned vertex UV preservation

Additionally, the existing `MultiLevelHierarchy` test infers that multiple levels were triggered
indirectly (via `setMinCoarseVertices(10)`) but never asserts the actual level count.

## New tests

### `HLSCMInternal.BuildHierarchy_LevelCount`
Call `buildHierarchy` directly on a small mesh with known vertex count and controlled
`minCoarseVerts` / `levelRatio`. Assert:
- `levels.size()` equals the expected count
- Each level's vertex count is approximately `prevCount / levelRatio`
- `level.localToOriginal` and `level.originalToLocal` are inverses of each other
- Pin vertices survive at every level

### `HLSCMInternal.ProlongateUVs_BarycentricReconstruction`
Construct a simple 2-level hierarchy manually (or via `buildHierarchy` on a 9-vertex grid).
Assign known UVs at the coarse level. Call `prolongateUVs`. Assert that each prolongated vertex's
UV matches the expected barycentric interpolation of its containingTri's UVs.

### `HLSCMInternal.SolveLSCMLevel_KnownMesh`
Call `solveLSCMLevel` directly on a pyramid `HierarchyLevel`. Assert:
- Returned UVs are finite and z=0
- The two pinned vertices have exactly the prescribed UV positions
- Result matches the top-level `AngleBasedLSCM::Compute` output on the same mesh

### `HLSCM.MultiLevelHierarchy` (update existing)
Add a way to observe the actual level count (e.g. expose it from `buildHierarchy` or via a
test hook). Replace the indirect `setMinCoarseVertices(10)` + validity check with an explicit
`EXPECT_GE(levelCount, 2)`.

## Acceptance criteria
- [ ] `HLSCMInternal.BuildHierarchy_LevelCount` passes and asserts level count directly
- [ ] `HLSCMInternal.ProlongateUVs_BarycentricReconstruction` passes
- [ ] `HLSCMInternal.SolveLSCMLevel_KnownMesh` passes and agrees with `AngleBasedLSCM`
- [ ] `HLSCM.MultiLevelHierarchy` directly asserts `levelCount >= 2`
- [ ] No changes to public API required (tests access `detail::hlscm` directly)
