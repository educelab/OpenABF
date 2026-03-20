# A2 — FacePtr range iteration requires awkward `*face` dereference

## GitHub Issue
https://github.com/educelab/OpenABF/issues/4

## Summary
Iterating over the edges of a face pointer requires `for (const auto& e : *face)`
rather than the more natural `for (const auto& e : face->edges())`. Adding a named
`edges()` accessor to `FacePtr` / `Face` improves usability significantly.

## File
`include/OpenABF/HalfEdgeMesh.hpp`

## Related
GitHub issue #4. This is subsumed by P1 — implement as part of P1.

## Acceptance Criteria
- [ ] `face->edges()` returns an iterable over the face's edges
- [ ] `for (const auto& e : face->edges())` works in range-based for loops
- [ ] `for (const auto& e : *face)` continues to work (no breaking change)
- [ ] Internal code in ABF/LSCM updated to use `face->edges()` where appropriate
- [ ] All existing tests pass

## Dependencies
- P1 (this is part of the P1 iteration redesign work)
