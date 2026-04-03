# A7 — Multi-pin UV constraints for LSCM

## GitHub Issue
https://github.com/educelab/OpenABF/issues/42

## Summary
`AngleBasedLSCM` supports pinning exactly two vertices. Allowing the caller to
supply three or more vertices with explicit UV positions enables boundary stitching
in multi-component pipelines (F3), fitting patches to an atlas, and reducing
distortion by adding anchors along the boundary.

## File
`include/OpenABF/AngleBasedLSCM.hpp`

## Acceptance Criteria
- [ ] `Compute(mesh, PinMap)` overload added where `PinMap` is
  `std::vector<std::pair<std::size_t, Vec<T,2>>>`
- [ ] `setPins(PinMap)` method added to the instance form
- [ ] All existing tests (default two-pin path and A5 explicit two-pin) pass unchanged
- [ ] Tests verify that pinned vertices land exactly at the specified UV positions

## Dependencies
- A5 (two-pin explicit override) must be merged first
