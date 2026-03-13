# M1 — `Face::barycenter()` hardcodes `Vec<T, 3>` instead of `Vec<T, Dim>`

## GitHub Issue
https://github.com/educelab/OpenABF/issues/25

## Summary
`Face::barycenter()` returns `Vec<T, 3>` regardless of the mesh's `Dim` template
parameter. It should use `Vec<T, Dim>`.

## File
`include/OpenABF/HalfEdgeMesh.hpp`

## Acceptance Criteria
- [ ] Return type changed to `Vec<T, Dim>`
- [ ] All tests pass

## Dependencies
None.
