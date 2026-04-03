# A6 — AngleBasedLSCM A/B matrix construction is repetitive and hard to follow

## GitHub Issue
https://github.com/educelab/OpenABF/issues/3

## Summary
The A and B matrix assembly in `AngleBasedLSCM::Compute` is a large block of nearly
identical per-vertex conditionals (repeated three times per face for e0, e1, e2).
Extracting a helper reduces duplication and makes the mathematical structure visible.

## File
`include/OpenABF/AngleBasedLSCM.hpp`

## Related
GitHub issue #3

## Acceptance Criteria
- [ ] The per-vertex triplet contribution is extracted into a helper lambda or inline
  function with a clear name
- [ ] The face loop body is visibly shorter and the structure matches the paper
- [ ] No behavioral change — all existing tests pass

## Dependencies
None.
