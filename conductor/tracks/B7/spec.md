# B7 — `PlanGrad` allocates and immediately discards the wheel vector

## GitHub Issue
https://github.com/educelab/OpenABF/issues/12

## Summary
In `ABF.hpp`, `PlanGrad` calls `v->wheel()` twice: once storing the result in an
unused `edges` variable, and again in the range-based for loop. This wastes a heap
allocation on every call.

## File
`include/OpenABF/ABF.hpp` (line 98)

## Acceptance Criteria
- [ ] The unused `auto edges = v->wheel();` line is removed
- [ ] All existing tests pass

## Dependencies
None. (Superseded in scope by P1, but this is a trivial isolated fix.)
