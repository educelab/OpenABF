# P5 — HLSCM hot-path allocation reduction

## GitHub Issue
https://github.com/educelab/OpenABF/issues/63

## Problem
Several data structures in the HLSCM decimation and solver pipeline use heap-heavy containers
where cache-friendly flat arrays would be faster. On meshes with 500K+ vertices these cause
millions of heap allocations per hierarchy level.

## Items (by priority)

### High
1. **UV map → dense vector** (`solveLSCMLevel`, `prolongateUVs`, output UV map):
   `unordered_map<size_t, array<T,2>>` → `vector<array<T,2>>` sized to original vertex count,
   with uninitialized entries for vertices not yet solved. Eliminates per-entry heap allocation
   and hash-table overhead; enables sequential access in prolongation.

2. **`vertexNeighbors()` post-collapse allocation:**
   Return the updated neighbor set as part of `CollapseResult` (or a separate out-param) so
   the PQ update after collapse reuses data already computed inside `tryCollapse`. Eliminates
   one O(valence) allocation per successful collapse.

### Medium
3. **`HierarchyLevel::originalToLocal` → dense vector:**
   Replace `unordered_map<size_t,size_t>` with `vector<size_t>` (SIZE_MAX sentinel) indexed by
   original vertex index. O(1) lookup with sequential allocation; size is known at level build time.

4. **`buildLevelMesh` face conversion:**
   Add an `insert_faces` overload accepting a `span`/range of `array<size_t,3>` to avoid
   constructing a `vector<vector<size_t>>` per face. Alternatively, accumulate into a flat buffer
   and call the existing overload once.

5. **`buildEdges_()` incremental maintenance:**
   Instead of rebuilding the full edge set from all faces at the start of each level, maintain
   edges incrementally during collapses using the neighbor results from `vertexNeighbors`. Remove
   the `buildEdges_()` call from `decimate()` loop body.

## Acceptance criteria
- [ ] No `unordered_map` in the UV storage path (`solveLSCMLevel`, `prolongateUVs`, output map)
- [ ] `vertexNeighbors()` does not perform a heap allocation when called immediately after `tryCollapse`
- [ ] `HierarchyLevel::originalToLocal` is a `vector<size_t>` with SIZE_MAX sentinel
- [ ] `buildLevelMesh` does not allocate a `vector<size_t>` per face
- [ ] `buildEdges_()` is not called inside the main decimation loop (or is eliminated entirely)
- [ ] All existing HLSCM tests pass
- [ ] No regression in parameterization quality on hemisphere/wavy/grid test meshes
