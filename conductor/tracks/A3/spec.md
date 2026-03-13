# A3 — `insert_face` boundary update behavior is undocumented for variadic form

## GitHub Issue
https://github.com/educelab/OpenABF/issues/16

## Summary
The variadic `insert_face(Args... args)` does not call `update_boundary()`, but this
is not clearly documented (unlike the vector form which has an explicit note). Users
expecting a single face insertion to be a complete operation may be surprised.

## File
`include/OpenABF/HalfEdgeMesh.hpp`

## Acceptance Criteria
- [ ] The variadic `insert_face` Doxygen comment explicitly states boundary is NOT
  updated and instructs callers to follow with `update_boundary()`
- [ ] The vector-form `insert_face(Vector&&)` comment is also consistent
- [ ] No behavioral changes; this is a documentation fix

## Dependencies
None.
