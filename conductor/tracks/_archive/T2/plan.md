# T2 Implementation Plan

## Phase 1 — Tests
In `tests/src/TestHalfEdgeMesh.cpp`, add `TEST(HalfEdgeMesh, FindPath_Disconnected)`:
1. Build a mesh with two disconnected triangles (vertices 0-2 in CC1, vertices 3-5 in CC2)
2. Call `FindEdgePath(mesh, 0, 3)` — source and destination in different CCs
3. Assert the returned vector is empty

## Phase 2 — Verify
- Run `ctest`
