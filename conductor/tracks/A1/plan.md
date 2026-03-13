# A1 Implementation Plan

## Phase 1 — Tests
1. Add a test that runs ABFPlusPlus with a tight threshold (e.g. 1e-6) and verifies
   the final gradient is ≤ threshold (or maxIters is hit)
2. Add a test that a loose threshold (e.g. 0.5) causes early termination

## Phase 2 — ABF
1. Add `T gradThreshold_{0.001}` member
2. Add `void setGradientThreshold(T t) { gradThreshold_ = t; }`
3. Update `compute()` to pass `gradThreshold_` to `Compute`
4. Update `static Compute(...)` signature to accept `gradThreshold = T(0.001)`
5. Replace hardcoded `0.001` with the parameter

## Phase 3 — ABFPlusPlus
Apply the same changes.

## Phase 4 — Verify
- Run `ctest`
- Confirm default behavior is unchanged
