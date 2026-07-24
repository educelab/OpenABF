# A10 — Add static Compute() overloads accepting levelRatio and minCoarseVertices

## GitHub Issue
https://github.com/educelab/OpenABF/issues/68

## Summary
Add static `Compute()` overloads to HLSCM that accept `levelRatio` and `minCoarseVertices` parameters directly, without requiring an instance to be constructed.

## Acceptance Criteria
- [ ] Static `Compute()` overloads accepting levelRatio and/or minCoarseVertices
- [ ] Consistent with existing static API patterns in the codebase
- [ ] Unit tests covering new overloads
