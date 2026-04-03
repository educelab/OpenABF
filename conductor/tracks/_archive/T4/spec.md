# T4 — No test for `Vec` reverse iterators (masks bug B3)

## GitHub Issue
https://github.com/educelab/OpenABF/issues/24

## Summary
The incorrect return types on `Vec::rbegin()`/`rend()` (B3) are not caught by the
test suite because no test exercises reverse iteration on `Vec`.

## Acceptance Criteria
- [ ] Tests exercise `rbegin()`/`rend()` and `crbegin()`/`crend()` on `Vec3f`
- [ ] Tests verify elements are visited in reverse order

## Dependencies
- Required by B3 (TDD: tests must be written before the fix)
