# T8 Implementation Plan

## Phase 1: Investigation
- [x] 1.1 Review InstanceAPILevelRatio test
- [x] 1.2 Identify observable effects of levelRatio parameter
  - `detail::hlscm::buildHierarchy` returns the `levels` vector directly;
    tests already access `detail::hlscm` (see `HLSCMInternal.DecimationMesh_RejectsPinnedVertex`).
    Calling it with two different `levelRatio` values on the same mesh gives a
    direct, observable measure (level count).

## Phase 2: Strengthen Test
- [x] 2.1 Add assertions that verify levelRatio affects hierarchy construction
  - Invoke `buildHierarchy` with `levelRatio=2` and `levelRatio=8` on a 20x20
    wavy mesh with `minCoarseVerts=10`, then `EXPECT_GT(small, large)`.
- [x] 2.2 Consider testing edge cases (very high/low ratios)
  - Used ratio=2 (minimum allowed, very fine) vs ratio=8 (aggressive) to keep
    the signal strong without overfitting; existing `MultiLevelHierarchy` test
    covers the small-mesh single-level path.

## Phase 3: Verification
- [x] 3.1 Run full test suite — all 6 ctest targets pass
