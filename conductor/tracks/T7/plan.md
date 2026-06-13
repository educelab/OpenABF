# T7 Implementation Plan

## Phase 1: Investigation
- [x] 1.1 Review ABFPlusPlusAnglePreservation test
- [x] 1.2 Determine if it passes vacuously on flat mesh
  - Wavy mesh (z = 0.3·sin·cos) produced UV maxDiff ≈ 0.057, L2 ≈ 0.62
    between ABF→HLSCM and geometry-only HLSCM. Non-zero, but weak — a regression
    that discarded ABF angles would only drop a small signal, and the existing
    `> 1e-6f` threshold could be tripped by unrelated numerical noise.
  - Test was not strictly vacuous, but the signal-to-noise ratio was poor.

## Phase 2: Fix
- [x] 2.1 Switch test to curved mesh
  - Replaced `ConstructWavySurface<…>(20, 20)` with
    `ConstructHemisphere<…>(12, 24)`. The hemisphere has genuine Gaussian
    curvature so ABF angle correction has real work to do.
- [x] 2.2 Tighten assertions
  - Snapshot raw geometry angles before `ABFType::Compute`, then assert
    `maxAngleDelta > 1e-3f` to prove ABF actually changed the input.
  - Assert UV `maxDiff > 1e-2f` (vs original `> 1e-6f`). On the hemisphere
    observed maxDiff ≈ 1.03 and L2 ≈ 9.08 — orders of magnitude above the
    threshold, but a regression that discards ABF angles drops both to ~0.

## Phase 3: Verification
- [x] 3.1 Run full test suite — all 6 ctest suites pass.
