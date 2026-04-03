# M4 — `operator<<` for `Vec` defined outside the `OpenABF` namespace

## GitHub Issue
https://github.com/educelab/OpenABF/issues/28

## Summary
The `operator<<` streaming operator for `Vec` is defined after the `OpenABF` namespace
closing brace. While ADL handles this in practice, it is cleaner to define it inside
the namespace or as an inline friend in the class body.

## File
`include/OpenABF/Vec.hpp` (line 274)

## Acceptance Criteria
- [ ] `operator<<` is moved inside the `OpenABF` namespace (or made a `friend` in `Vec`)
- [ ] All tests pass
- [ ] No change in observable behavior

## Dependencies
None.
