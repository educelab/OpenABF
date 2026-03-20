# B2 — ABFPlusPlus LambdaStarInv inversion uses hardcoded `1.F` (float)

## GitHub Issue
https://github.com/educelab/OpenABF/issues/7

## Summary
In `ABFPlusPlus::Compute`, the LambdaStarInv matrix element inversion uses the
float literal `1.F` regardless of the template scalar type `T`. When `T = double`,
results are silently truncated to single precision.

## File
`include/OpenABF/ABFPlusPlus.hpp` (line 187)

## Acceptance Criteria
- [ ] `1.F / it.value()` replaced with `T(1) / it.value()`
- [ ] All existing tests pass
- [ ] A `double`-precision ABFPlusPlus test does not lose accuracy (see F4)

## Dependencies
None.
