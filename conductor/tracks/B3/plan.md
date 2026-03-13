# B3 Implementation Plan

## Phase 1 — Tests
1. In `tests/src/TestVec.cpp`, add a `TEST(Vec, ReverseIteration)` case that:
   - Constructs a `Vec3f`
   - Iterates with `rbegin()`/`rend()` and verifies element order is reversed
   - Uses `crbegin()`/`crend()` on a const `Vec3f`

## Phase 2 — Fix
1. Change `rbegin()` return type from `iterator` to `reverse_iterator`
2. Change `rbegin() const` return type from `const_iterator` to `const_reverse_iterator`
3. Change `crbegin()` return type from `const_iterator` to `const_reverse_iterator`
4. Repeat for `rend()`, `rend() const`, `crend()`

## Phase 3 — Verify
- Compile and run `ctest`
- Confirm new reverse iteration test passes
