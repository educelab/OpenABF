# B2 Implementation Plan

## Phase 1 — Fix
1. In `ABFPlusPlus::Compute`, locate the LambdaStarInv inversion loop
2. Replace `it.valueRef() = 1.F / it.value();` with `it.valueRef() = T(1) / it.value();`

## Phase 2 — Verify
- Run `ctest`
- Optionally add a double-precision smoke test (or defer to F4 track)
