# P5 Implementation Plan

## Phase 1: UV map → dense vector
- [ ] 1.1 Change `solveLSCMLevel` return type from `unordered_map<size_t, array<T,2>>` to `vector<array<T,2>>` (size = original vertex count, sentinel = `{NaN, NaN}`)
- [ ] 1.2 Update `prolongateUVs` to accept/return `vector<array<T,2>>`
- [ ] 1.3 Update `buildInitialGuess` and UV extraction loop to use the vector
- [ ] 1.4 Update `ComputeImpl` UV output loop
- [ ] 1.5 Run tests — verify no regression

## Phase 2: vertexNeighbors post-collapse
- [ ] 2.1 Modify `tryCollapse` to return neighbors of `vKeep` as part of its result (or a separate overload)
- [ ] 2.2 Update the collapse loop in `buildHierarchy` to pass the returned neighbors to PQ update instead of calling `vertexNeighbors()`
- [ ] 2.3 Remove or make private the standalone `vertexNeighbors()` method if no longer needed externally

## Phase 3: HierarchyLevel::originalToLocal
- [ ] 3.1 Change `HierarchyLevel::originalToLocal` from `unordered_map<size_t,size_t>` to `vector<size_t>` with `SIZE_MAX` sentinel
- [ ] 3.2 Update `snapshot()` to fill the vector (size = finest-mesh vertex count, passed in as param)
- [ ] 3.3 Update all lookup sites in `solveLSCMLevel`

## Phase 4: buildLevelMesh face conversion
- [ ] 4.1 Add `insert_faces(span<const array<size_t,3>>)` overload to `HalfEdgeMesh` (or use a flat `vector<array<size_t,3>>` buffer passed as a single call)
- [ ] 4.2 Update `buildLevelMesh` to use the new overload

## Phase 5: buildEdges_ incremental (optional, largest scope)
- [ ] 5.1 Maintain `edges_` set incrementally in `executeCollapse` using neighbor data
- [ ] 5.2 Remove `buildEdges_()` call from level initialization

## Phase 6: Verification
- [ ] 6.1 `cmake --build` clean build
- [ ] 6.2 `ctest` all tests pass
- [ ] 6.3 Profile on a 500K-vertex mesh before/after to confirm allocation reduction
