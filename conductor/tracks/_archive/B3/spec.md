# B3 — `Vec` reverse iterators have wrong return type

## GitHub Issue
https://github.com/educelab/OpenABF/issues/8

## Summary
`Vec::rbegin()`, `rend()`, `crbegin()`, `crend()` are declared to return `iterator`
/ `const_iterator` but must return `reverse_iterator` / `const_reverse_iterator`.
This causes silent type mismatches and prevents correct reverse iteration.

## File
`include/OpenABF/Vec.hpp` (lines 110–127)

## Acceptance Criteria
- [ ] Tests exist that exercise reverse iteration on `Vec` (see T4)
- [ ] All four reverse iterator methods return the correct iterator type
- [ ] Range-based reverse iteration (`for (auto it = v.rbegin(); it != v.rend(); ++it)`)
  compiles and produces correct results

## Dependencies
- T4 should be added first (TDD)
