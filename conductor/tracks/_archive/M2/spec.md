# M2 — `Edge::magnitude()` is non-const

## GitHub Issue
https://github.com/educelab/OpenABF/issues/26

## Summary
`Edge::magnitude()` only reads vertex positions but is not declared `const`, preventing
its use on const-qualified edge references. `Face::area()` (which is const) depends on
it via a non-const path through a shared_ptr, which obscures the issue.

## File
`include/OpenABF/HalfEdgeMesh.hpp`

## Acceptance Criteria
- [ ] `magnitude()` declared `const`
- [ ] All tests pass

## Dependencies
None.
