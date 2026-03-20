# F4 Implementation Plan

## Phase 1 — Tests
1. Add `TEST(Parameterizations, ABF_double)`, `TEST(Parameterizations, ABFPlusPlus_double)`,
   and `TEST(Parameterization, AngleBasedLSCM_double)` using the same pyramid fixture
   but with `double` scalar type
2. Assert double results are equal-or-better than float (tighter gradient, more precise UV)

## Phase 2 — Verify
- Run `ctest`
- Confirm B2 is resolved first; double tests will fail otherwise
