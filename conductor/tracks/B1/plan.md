# B1 Implementation Plan

## Phase 1 — Tests
1. Add a flat grid mesh fixture to the test suite (see T1 for shared work)
2. Add `TEST(Parameterizations, ABF_Grid)` that verifies correct UV output for a mesh
   with ≥ 3 interior vertices using `ABF::Compute` + `AngleBasedLSCM::Compute`

## Phase 2 — Fix
1. In `ABF::Compute`, before the interior-vertex lambda update loop, capture a fixed
   base index: `const auto lambdaBase = edgeCnt + faceCnt;`
2. Replace `delta(idx + intIdx, 0)` → `delta(lambdaBase + intIdx, 0)`
3. Replace `delta(idx + vIntCnt + intIdx, 0)` → `delta(lambdaBase + vIntCnt + intIdx, 0)`
4. Remove the `idx++` inside the loop (it is no longer needed)

## Phase 3 — Verify
- Run full test suite: `ctest`
- Confirm new grid test passes
- Confirm all existing parameterization tests pass
