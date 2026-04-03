# T7 — ABFPlusPlusAnglePreservation may be vacuous on flat wavy mesh

## GitHub Issue
https://github.com/educelab/OpenABF/issues/50

## Summary
The ABFPlusPlusAnglePreservation test may be vacuously passing on the flat wavy mesh because ABF angle preservation is not meaningfully tested by a flat mesh. Strengthen or replace with a curved mesh that actually exercises ABF angle constraints.

## Acceptance Criteria
- [ ] Investigate whether current test is vacuous
- [ ] Replace or supplement with a curved mesh test
- [ ] Test meaningfully verifies ABF angle preservation
