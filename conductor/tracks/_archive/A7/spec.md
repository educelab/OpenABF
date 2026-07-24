# A7 — Multi-pin UV constraints for LSCM and Hierarchical LSCM

## GitHub Issue
https://github.com/educelab/OpenABF/issues/42

## Summary
`AngleBasedLSCM` and `HierarchicalLSCM` currently support pinning exactly two
vertices. Allowing the caller to supply three or more vertices with explicit UV
positions enables boundary stitching in multi-component pipelines (F3), fitting
patches to an atlas, and reducing distortion by adding anchors along the
boundary.

Multi-pin support is added to **both** solvers via a shared assembly helper.
Pins are guaranteed to survive every decimation level in HLSCM (the existing
`isPinned_` mechanism already protects them from collapse), so HLSCM multi-pin
is achievable without algorithmic changes — only API surface and per-level
bookkeeping.

## Files
- `include/OpenABF/detail/LSCMSystem.hpp` — new `buildSystemMultiPin(mesh, PinMap)`
  overload that accepts caller-supplied UVs (no axis-snap auto-placement).
- `include/OpenABF/AngleBasedLSCM.hpp` — multi-pin `Compute` overload + `setPins`.
- `include/OpenABF/HierarchicalLSCM.hpp` — multi-pin `Compute` overload + `setPins`;
  generalize `DecimationMesh::build` to flag all pins as non-collapsible;
  generalize `solveLSCMLevel` to map a PinMap to level-local indices.

## Acceptance Criteria
- [ ] `using PinMap = std::vector<std::pair<std::size_t, Vec<T,2>>>` defined and
  re-exported by both `AngleBasedLSCM` and `HierarchicalLSCM`
- [ ] `AngleBasedLSCM::Compute(mesh, PinMap)` static overload + `setPins(PinMap)`
  instance method
- [ ] `HierarchicalLSCM::Compute(mesh, PinMap)` static overload + `setPins(PinMap)`
  instance method
- [ ] All existing tests (default two-pin path and A5 explicit two-pin) pass
  unchanged for both solvers
- [ ] Tests verify that pinned vertices land exactly at the specified UV
  positions for both solvers (≥3 pins)
- [ ] Tests cover the instance form (`setPins() + compute()`) for both solvers
- [ ] HLSCM with many pins still builds a valid hierarchy (or gracefully falls
  back to single-level LSCM when decimation is exhausted)

## Dependencies
- A5 (two-pin explicit override) merged
- F1 (HLSCM) merged
