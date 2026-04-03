# T8 — Strengthen InstanceAPILevelRatio to verify parameters are applied

## GitHub Issue
https://github.com/educelab/OpenABF/issues/49

## Summary
The InstanceAPILevelRatio test verifies that the API call succeeds but does not verify that the levelRatio parameter actually affects the hierarchy construction. Strengthen it to verify observable effects.

## Acceptance Criteria
- [ ] Test verifies that different levelRatio values produce different hierarchy depths or coarse vertex counts
- [ ] Test is not vacuous — a broken levelRatio parameter would cause it to fail
