# P5 Implementation Plan

Branch: `p5-hlscm-alloc-reduction` (in `../OpenABF-p5` worktree).

## Phase 1: UV map → dense vector — committed
- [x] 1.1 Change `solveLSCMLevel` return type to `vector<optional<Vec<T,2>>>` (aliased as `UVVector<T>`) sized to the finest-mesh vertex count; unset slots hold `std::nullopt`.
- [x] 1.2 Update `prolongateUVs` signature accordingly (input vector is moved in, output vector is moved out). Body uses `Vec<T,2>` arithmetic for the barycentric mix.
- [x] 1.3 Update `buildInitialGuess` to dispatch on `optional::has_value()`.
- [x] 1.4 Update `ComputeImpl` final UV output loop — every finest-level vertex has a UV so the unwrap is safe.
- [x] 1.5 ctest — all 43 parameterization tests pass.

## Phase 2: vertexNeighbors post-collapse — committed
- [x] 2.1 `tryCollapse` accepts optional `std::vector<std::size_t>* outKeepNbrs`; populates it with vKeep's post-collapse neighbors using already-compacted `vertFaces_[vKeep]`.
- [x] 2.2 `buildHierarchy` reuses a single `nbrs` vector across collapses, passed by pointer to `tryCollapse`; eliminates one heap allocation per successful collapse.
- [x] 2.3 `vertexNeighbors(v)` retained as a public read-only accessor — still used by the test suite and harmless when not on the hot path.

## Phase 3: HierarchyLevel::originalToLocal — committed
- [x] 3.1 Changed `HierarchyLevel::originalToLocal` to `vector<optional<size_t>>`; absent vertices hold `std::nullopt` (type system enforces unwrap rather than reserving a sentinel `size_t` value).
- [x] 3.2 `snapshot()` sizes the vector to `alive_.size()` (finest-mesh vertex count) on entry and `assign`s `std::nullopt`.
- [x] 3.3 Lookups in `solveLSCMLevel` use `*operator[]`; tests check `has_value()`.

## Phase 4: buildLevelMesh face conversion — committed
- [x] 4.1 No new overload needed — `HalfEdgeMesh::insert_faces` is already generic over containers-of-iterables and accepts `vector<array<size_t,3>>` directly.
- [x] 4.2 `buildLevelMesh` now passes `level.faces` straight to `insert_faces`, eliminating the per-face `vector<vector<size_t>>` heap allocation.

## Phase 5: buildEdges_ incremental — TRIED AND DROPPED
Implemented in commit `886f14c` then reverted before review. The implementation maintained `edges_` incrementally inside `tryCollapse` (unordered_set of packed `min(a,b)*N + max(a,b)` keys, with a lazy vector view for the PQ-seed loop) and removed the per-level `buildEdges_` rebuild. Tests passed (44/44 parameterization assertions, including a new invariant test).

Profiled on the wavy built-in series 50k–800k against the phase-1-4 binary (`baf00bb`): HLSCM-LSCG and HLSCM-CG deltas were all within ±2.5%, indistinguishable from run-to-run noise.

Why: `buildEdges_()` runs ~5–7 times per `Compute()` (once per hierarchy level), totalling ~10 ms even on a 1M-face mesh. HLSCM-CG at 1M spends ~459 s in the CG solves at each level. Eliminating buildEdges_ is at most a 0.003% improvement on the overall HLSCM solve — well below the noise floor, while adding a non-trivial invariant for `tryCollapse` to maintain. Not worth the maintenance cost.

If a future profiling pass shows `buildEdges_` becoming a meaningful share of the HLSCM budget (e.g. after the CG solve gets dramatically faster from a different optimization), revisit the experiment from `886f14c` in the git history.

## Phase 6: Verification — done
- [x] 6.1 Clean ninja build of P5 worktree (`build/`).
- [x] 6.2 ctest — 6/6 test binaries pass, 43/43 parameterization assertions pass.
- [x] 6.3 Profiled phases 1–4 vs the post-A8 baseline on the wavy built-in series 50k–1M at 1 thread. Deltas in the HLSCM columns are within run-to-run noise at every size — the alloc-pressure savings phases 1–4 deliver don't show up against the dominant CG-solve cost on these meshes. Allocator wins should be more visible under memory-pressure workloads (e.g. running many parallel HLSCM Compute()s) or after the per-level CG solve gets faster.

## Phase 7: Conductor & GitHub
- [x] 7.1 Commits pushed (phases 1–4 only): `a14d34a`, `baf00bb`.
- [x] 7.2 PR #92 open against #63.
- [x] 7.3 Update `tracks.md` (move P5 to archived once PR lands).
