# T2 — No test for `FindEdgePath` on disconnected mesh

## GitHub Issue
https://github.com/educelab/OpenABF/issues/22

## Summary
`FindEdgePath` returns an empty vector when no path exists between two vertices
(including when they are in different connected components), but this case is untested.

## Acceptance Criteria
- [ ] A test calls `FindEdgePath` on a mesh with two disconnected components where
  source and destination are in different components
- [ ] The test asserts the returned path is empty

## Dependencies
None.
