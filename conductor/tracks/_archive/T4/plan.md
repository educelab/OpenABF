# T4 Implementation Plan

## Phase 1 — Tests
In `tests/src/TestVec.cpp`, add `TEST(Vec, ReverseIteration)`:
1. Construct `Vec3f v{1.f, 2.f, 3.f}`
2. Collect elements using `rbegin()`/`rend()` into a vector
3. Assert the collected order is `{3.f, 2.f, 1.f}`
4. Repeat for `crbegin()`/`crend()` on a const Vec3f

## Phase 2 — Verify
- Confirm test fails to compile before B3 fix (return type mismatch)
- After B3 fix, run `ctest` and confirm test passes
