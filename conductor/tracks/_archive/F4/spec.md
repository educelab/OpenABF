# F4 — Double-precision test coverage for parameterization algorithms

## GitHub Issue
https://github.com/educelab/OpenABF/issues/20

## Summary
All parameterization tests use `float`. Given the library's stated goal of high
numerical accuracy for research use, tests verifying that `double` produces tighter
results than `float` on the same mesh would strengthen confidence.

## Acceptance Criteria
- [ ] `ABF`, `ABFPlusPlus`, and `AngleBasedLSCM` tests include `double`-precision
  variants alongside the existing `float` variants
- [ ] Double-precision results are verified to be at least as accurate as float
  (e.g. lower final gradient, closer to analytical solution)
- [ ] B2 (the `1.F` precision bug) must be fixed first

## Dependencies
- B2 must be fixed before double-precision results can be trusted
