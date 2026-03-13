# T1 — ABF/ABFPlusPlus only tested on a single-interior-vertex pyramid

## GitHub Issue
https://github.com/educelab/OpenABF/issues/21

## Summary
The pyramid mesh used in all parameterization tests has exactly one interior vertex,
which masks any bug that requires 2+ interior vertices (including B1). Tests must
cover meshes with more interior vertices and edge cases.

## Acceptance Criteria
- [ ] A shared flat-grid mesh fixture (e.g. 3×3 grid = 4 interior vertices) is added
  to `tests/src/Utils.hpp`
- [ ] `ABF`, `ABFPlusPlus`, and `AngleBasedLSCM` are each tested on this grid mesh
- [ ] A test for a mesh with 0 interior vertices (single triangle or boundary-only)
  does not crash
- [ ] A test verifies behavior when `maxIters` is reached without convergence

## Dependencies
- Required by B1 (TDD: tests must be written before the fix)
