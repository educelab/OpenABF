# P6 — Precompute sin/cos values before HLSCM face assembly loop

## GitHub Issue
https://github.com/educelab/OpenABF/issues/66

## Summary
Precompute trigonometric sin/cos values used in the HLSCM face assembly loop to avoid redundant computation per iteration.

## Acceptance Criteria
- [ ] Identify all sin/cos calls inside HLSCM face assembly loop
- [ ] Precompute values before the loop and store in local arrays
- [ ] Verify no regression in parameterization correctness
- [ ] Benchmark to confirm performance improvement
