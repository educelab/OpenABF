# F3 Implementation Plan

Status: completed (2026-06-12)

The merged design replaces the originally-planned free function
`ExtractConnectedComponents` in `HalfEdgeMeshUtils.hpp` with a method on
`HalfEdgeMesh` that returns an `ExtractedComponent` struct (with
`vertex_map` and `face_map`). The high-level `ParameterizeConnectedComponents`
wrapper (original Phase 3) was descoped — see note below.

## Phase 1 — Tests
- [x] **Task 1.1**: Test `extract_connected_components` on a single-CC mesh
  (returns one element, same geometry) `a32e1e2`
- [x] **Task 1.2**: Test on a two-CC mesh (returns two elements, correct
  vertex maps) `a32e1e2`
- [x] **Task 1.3**: Test trait preservation on a mesh with pre-computed
  angle traits `a32e1e2`

## Phase 2 — extract_connected_components
- [x] **Task 2.1**: Implement `HalfEdgeMesh::extract_connected_components()`
  as a method (not a free function in `HalfEdgeMeshUtils.hpp`) `a32e1e2`
- [x] **Task 2.2**: Return `ExtractedComponent` struct with `mesh`,
  `vertex_map`, and `face_map` fields `c97395f`

## Phase 3 — ParameterizeConnectedComponents
- [-] **Task 3.1**: Convenience wrapper — **DESCOPED** (2026-06-12). A
  bundled wrapper that fixes the angle-optimizer + parameterizer choice and
  writes UVs back automatically doesn't carry its weight: callers already
  have the `vertex_map` they need, and per-component pipeline choices
  (which optimizer, which solver, pin selection, error handling) vary too
  much to hide behind one helper. The low-level extractor is sufficient.

## Phase 4 — Verify
- [x] **Verify 4.1**: `ctest` green on develop after merges `a32e1e2`, `c97395f`
- [x] **Verify 4.2**: Reviewed and merged into develop via PRs #80, #81
