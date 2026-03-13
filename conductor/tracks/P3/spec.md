# P3 — `InitializeAnglesAndWeights` calls `v->wheel()` redundantly

## GitHub Issue
https://github.com/educelab/OpenABF/issues/14

## Summary
In `InitializeAnglesAndWeights`, the interior vertex loop stores the result of
`v->wheel()` in a local variable and then passes over it with `std::accumulate`
followed by a second range-based for loop. This is two passes over the same data,
which is fine, but the variable capture pattern could be simplified.

## File
`include/OpenABF/ABF.hpp` (line 73)

## Acceptance Criteria
- [ ] `v->wheel()` is called at most once per interior vertex in this function
- [ ] All existing tests pass

## Dependencies
- P1: once wheel() returns a lazy view rather than a vector, this becomes a non-issue;
  P3 may be deferred or closed by P1
