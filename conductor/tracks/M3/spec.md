# M3 — `detail::erase_if` is dead code

## GitHub Issue
https://github.com/educelab/OpenABF/issues/27

## Summary
`detail::erase_if` in `HalfEdgeMesh.hpp` takes its container argument by value,
modifies the copy, and returns it. No call sites exist in the project — it is unused
dead code.

## File
`include/OpenABF/HalfEdgeMesh.hpp`

## Acceptance Criteria
- [ ] The function is removed
- [ ] All tests pass

## Dependencies
None.
