# A1 — Convergence tolerance is hardcoded in ABF and ABFPlusPlus

## GitHub Issue
https://github.com/educelab/OpenABF/issues/15

## Summary
The gradient convergence thresholds (`gradient > 0.001` and `gradDelta > 0.001`) are
hardcoded in both `ABF::Compute` and `ABFPlusPlus::Compute`. Users cannot tune
precision/speed trade-offs for their specific mesh sizes and use cases.

## Files
`include/OpenABF/ABF.hpp`, `include/OpenABF/ABFPlusPlus.hpp`

## Acceptance Criteria
- [ ] A `setGradientThreshold(T)` method is added alongside `setMaxIterations`
- [ ] The static `Compute` overloads accept an optional `gradThreshold` parameter
  (with default 0.001 for backward compatibility)
- [ ] Tests verify that tightening the threshold produces a smaller final gradient
- [ ] All existing tests pass with the default threshold

## Dependencies
None.
