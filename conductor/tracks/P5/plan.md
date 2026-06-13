# P5 Implementation Plan

Branch: `p5-hlscm-alloc-reduction` (in `../OpenABF-p5` worktree).

## Phase 1: UV map → dense vector — committed
- [x] 1.1 Change `solveLSCMLevel` return type to `vector<array<T,2>>` sized to the finest-mesh vertex count with `kUnsetUV<T>` (NaN) sentinel.
- [x] 1.2 Update `prolongateUVs` signature accordingly (input vector is moved in, output vector is moved out).
- [x] 1.3 Update `buildInitialGuess` to dispatch on `std::isnan(uv[0])` instead of `unordered_map::find`.
- [x] 1.4 Update `ComputeImpl` final UV output loop — every finest-level vertex has a UV so no `find()` is needed.
- [x] 1.5 ctest — all 43 parameterization tests pass.

## Phase 2: vertexNeighbors post-collapse — on disk, uncommitted
- [x] 2.1 `tryCollapse` accepts optional `std::vector<std::size_t>* outKeepNbrs`; populates it with vKeep's post-collapse neighbors using already-compacted `vertFaces_[vKeep]`.
- [x] 2.2 `buildHierarchy` reuses a single `nbrs` vector across collapses, passed by pointer to `tryCollapse`; eliminates one heap allocation per successful collapse.
- [x] 2.3 `vertexNeighbors(v)` retained as a public read-only accessor — still used by the test suite and harmless when not on the hot path.

## Phase 3: HierarchyLevel::originalToLocal — on disk, uncommitted
- [x] 3.1 Changed `HierarchyLevel::originalToLocal` to `vector<size_t>` with `HierarchyLevel::kAbsent = SIZE_MAX` sentinel.
- [x] 3.2 `snapshot()` sizes the vector to `alive_.size()` (finest-mesh vertex count) on entry.
- [x] 3.3 Lookups in `solveLSCMLevel` use `operator[]`; tests updated to use the sentinel-aware `isPresent()` helper.

## Phase 4: buildLevelMesh face conversion — on disk, uncommitted
- [x] 4.1 No new overload needed — `HalfEdgeMesh::insert_faces` is already generic over containers-of-iterables and accepts `vector<array<size_t,3>>` directly.
- [x] 4.2 `buildLevelMesh` now passes `level.faces` straight to `insert_faces`, eliminating the per-face `vector<vector<size_t>>` heap allocation.

## Phase 5: buildEdges_ incremental — DEFERRED
Inspection of the call sites shows `buildEdges_()` runs only **once per hierarchy level** (5–7 calls total on a 1M-face mesh), not per collapse. It's not on the hot path; the spec flagged this phase as "optional, largest scope" and the ROI is small relative to the disruption. Re-evaluate after profiling on a 500k+ mesh if needed.

## Phase 6: Verification — done
- [x] 6.1 Clean ninja build of P5 worktree (`build/`).
- [x] 6.2 ctest — 6/6 test binaries pass, 43/43 parameterization assertions pass.
- [x] 6.3 Smoke timing against P4-branch (post-A8) baseline at 50k/100k wavy: HLSCM-CG drops ~2–3 % at this size (within noise). Phase 1–4 allocator wins are designed to scale with mesh size; expect more on 500k+ meshes — recommend re-profiling once the user runs against real scroll data.

## Phase 7: Conductor & GitHub — pending commits/signing
- [ ] 7.1 Sign and push the P5 commits (phase 1 already committed; phases 2–4 + amalgamation staged).
- [ ] 7.2 Open PR against #63 referencing phases done and the phase-5 deferral.
- [ ] 7.3 Update tracks.md (move P5 to archived once PR lands).
