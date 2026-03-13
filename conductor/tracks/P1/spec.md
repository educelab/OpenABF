# P1 — Mesh iteration methods return `std::vector` by value (excessive allocation)

## GitHub Issue
https://github.com/educelab/OpenABF/issues/4

## Summary
Every mesh accessor — `vertices()`, `edges()`, `faces()`, `vertices_interior()`,
`vertices_boundary()`, `wheel()`, `boundaries()`, `outgoing_edges()`,
`connected_components()` — allocates and returns a `std::vector` by value on every
call. The ABF/ABFPlusPlus solvers call many of these in a tight per-iteration loop,
causing O(n) heap allocations per solver step that are immediately discarded.

## File
`include/OpenABF/HalfEdgeMesh.hpp`

## Related
GitHub issue #4

## Acceptance Criteria
- [ ] Hot-path iteration (used inside ABF/LSCM solve loops) does not allocate a new
  vector on every call
- [ ] Range-based for loops over vertices, edges, faces, and wheel remain syntactically
  identical or cleaner for callers
- [ ] A2 (FacePtr iteration) is resolved as part of this work
- [ ] All existing tests pass
- [ ] Solver performance on a large mesh is measurably improved (benchmark or profile)

## Proposed Approach
Replace vector-returning methods with lazy range views. Options:
1. **Span/view over internal storage** for `vertices()`, `edges()`, `faces()` (simplest
   since the internal `verts_` and `faces_` vectors are already contiguous)
2. **Custom iterator/sentinel** for `wheel()`, `outgoing_edges()`, `boundaries()` which
   are inherently linked-list traversals
3. **C++20 `std::ranges`** if the standard version can be bumped; otherwise a lightweight
   custom range type

A `face->edges()` named method returning an iterable (resolving A2) should be added
as part of this work.

## Dependencies
- A2 is subsumed by this track
- P2 and P3 may become simpler or redundant once this is done
- F3 benefits from this being completed first
