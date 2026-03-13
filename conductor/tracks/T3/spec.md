# T3 — No round-trip IO test for OBJ and PLY formats

## GitHub Issue
https://github.com/educelab/OpenABF/issues/23

## Summary
Reading a mesh, writing it, and reading it back should produce an identical mesh.
No such round-trip test exists for either OBJ or PLY.

## Acceptance Criteria
- [ ] OBJ round-trip: write a mesh to a string stream, read it back, assert vertex
  positions and face connectivity are identical
- [ ] PLY round-trip: same
- [ ] Tests use in-memory streams (no temp files required)

## Dependencies
None.
