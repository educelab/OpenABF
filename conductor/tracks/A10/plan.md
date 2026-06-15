# A10 Implementation Plan

## Phase 1: API Design
- [x] 1.1 Review existing static Compute() signatures
- [x] 1.2 Design new overload signatures

Note: The original issue (#68) proposed two new overloads:

  ```cpp
  Compute(mesh, levelRatio, minCoarseVerts);
  Compute(mesh, pin0, pin1, levelRatio, minCoarseVerts);
  ```

After PR #93 introduced the `PinMap` interface, the second shape was
adapted to `Compute(mesh, const PinMap&, levelRatio, minCoarseVerts)`.

The first shape (`Compute(mesh, levelRatio, minCoarseVerts)`) cannot be
added today because its `(Mesh::Pointer&, size_t, size_t)` signature
collides with the still-present `[[deprecated]] Compute(mesh, pin0Idx,
pin1Idx)`. It is blocked on the deprecated-overload cleanup tracked in
a separate issue (see #68 "blocked by" link). Until then, auto-pin
users who need custom tuning use the instance API (`set_level_ratio`,
`set_min_coarse_vertices`).

## Phase 2: Implementation
- [x] 2.1 Implement static overload in HLSCM (PinMap variant)
- [x] 2.2 Update single-header via amalgamation script

## Phase 3: Testing
- [x] 3.1 Write unit tests for new overload
- [x] 3.2 Run full test suite — all tests pass
