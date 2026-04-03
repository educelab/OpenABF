# A4 — `gradient()` not `[[nodiscard]]`; inconsistent with `iterations()`

## GitHub Issue
https://github.com/educelab/OpenABF/issues/17

## Summary
`ABF::iterations()` and `ABFPlusPlus::iterations()` are marked `[[nodiscard]]` but
`gradient()` is not, despite both returning state that is only meaningful after
`compute()`. Applying `[[nodiscard]]` consistently prevents silent misuse.

## Files
`include/OpenABF/ABF.hpp`, `include/OpenABF/ABFPlusPlus.hpp`

## Acceptance Criteria
- [ ] `gradient()` is marked `[[nodiscard]]` in both `ABF` and `ABFPlusPlus`
- [ ] All tests pass

## Dependencies
None.
