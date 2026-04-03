# T6 — Add test for HLSCM::setPinnedVertices() instance API

## GitHub Issue
https://github.com/educelab/OpenABF/issues/51

## Summary
Add a test that exercises the `HLSCM::setPinnedVertices()` instance API method to ensure it correctly constrains the UV parameterization.

## Acceptance Criteria
- [ ] Test creates a mesh with known topology
- [ ] Test calls setPinnedVertices() before Compute()
- [ ] Test verifies pinned vertices are fixed in output UV map
