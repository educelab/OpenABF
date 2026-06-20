# B8 — LSCM area-preserving rescale of flattening output

## GitHub Issue
https://github.com/educelab/OpenABF/issues/98

## Type
Bug (quality / numerical-accuracy)

## Summary
LSCM fixes the global scale of its output entirely from the two auto-selected
pins: `detail::lscm::AutoPlacePair` places pin1 at signed distance `|p1 - p0|`,
the raw **3D length of a single boundary edge**. Because LSCM is conformal — it
recovers the parameterization only up to a global similarity — that one edge
length is the sole determinant of the result's overall scale. If the
auto-selected boundary edge is anomalously stretched or compressed relative to
the rest of the surface, the entire flattening is globally mis-scaled, working
against the goal of preserving the original 3D surface area.

This track adds an **opt-in, post-solve uniform rescale** of the computed UVs
that better preserves the original surface area, shared by both
`AngleBasedLSCM` and `HierarchicalLSCM`.

## Problem Description
- `AutoPlacePair` (`include/OpenABF/detail/LSCMSystem.hpp`) bakes the raw 3D
  pin-edge length into pin1's UV, setting the global scale from one edge.
- `AutoSelectPins` picks the *first* boundary vertex and a boundary-adjacent
  neighbor with no regard for how representative that edge's length is.
- Both `AngleBasedLSCM::ComputeImpl` and `HierarchicalLSCM` consume this pin
  path and overwrite each vertex's 3D position with its UV, so both inherit the
  mis-scale.
- Distinct from the closed A5 (#2), which addressed *which* vertices to pin and
  *where*; this track corrects the **global scale of the result after
  flattening**, independent of pin choice.

## Design

### Rescale modes (selectable)
A new mode enum, defined in `detail::lscm` and re-exported as a public alias by
both solver classes:

- `Off` — no rescale. Current behavior; the right choice when the caller
  supplies explicit pins that already define the intended scale/placement.
- `LowStretchReference` — select the single least-distorted face (the most
  locally isometric: Jacobian singular values closest to equal / minimal shape
  distortion) and scale all UVs so that face's UV area equals its 3D area. Close
  to the issue's original phrasing; simpler but anchored on one face.
- `GlobalAreaMin` — choose one uniform factor that minimizes total area
  distortion across the whole mesh, robust to any single anomalous edge/face.
  Implementation sets the **area-weighted median per-face area ratio**
  (`a_3d / a_uv`) to 1, so a handful of degenerate faces cannot dominate the
  scale. (A total-area-matching factor `s = sqrt(Σa_3d / Σa_uv)` is a simpler
  alternative considered; the median form is preferred for outlier robustness.)

A uniform scale is applied about the UV origin to every vertex (`x`, `y`).

### API
- `set_rescale_mode(RescaleMode)` instance setter (+ private member, applied by
  `compute()`), on both `AngleBasedLSCM` and `HierarchicalLSCM`.
- Static `Compute` overloads gain an optional trailing
  `RescaleMode mode = RescaleMode::Off` parameter (both the auto-pin and
  explicit-`PinMap` forms).
- New public type alias `using RescaleMode = detail::lscm::RescaleMode;` on both
  classes (mirrors the existing `PinMap` re-export).

### Default and v3.0
- Default is `Off` for the current (2.x) line so existing output is unchanged.
- Doxygen must document that the **default changes to `GlobalAreaMin` in 3.0**.
  Coordinate the default flip with the 3.0 deprecation cleanup (#96).

### Implementation note (ordering)
The rescale needs the original 3D geometry, but `ComputeImpl` overwrites each
vertex's `pos` with its UV. The helper must therefore **snapshot original 3D
face areas before the UV writeback**, then compute UV areas from the written
UVs and apply the factor. For HLSCM the rescale is applied **once** to the final
top-level result.

### Caveat
A uniform scale about the origin moves pinned vertices off their specified UVs,
so rescaling is generally only meaningful on the auto-pin path. This must be
documented; `Off` remains the default precisely so explicit-pin callers are
unaffected.

## Files
- `include/OpenABF/detail/LSCMSystem.hpp` — `RescaleMode` enum; original-area
  snapshot helper; `RescaleUVs(mesh, mode, originalAreas)` shared helper.
- `include/OpenABF/AngleBasedLSCM.hpp` — `RescaleMode` alias, `set_rescale_mode`,
  static overload mode params, thread mode through `ComputeImpl`.
- `include/OpenABF/HierarchicalLSCM.hpp` — same surface; apply rescale once to
  the final result.
- `single_include/OpenABF/OpenABF.hpp` — regenerate via amalgamation script.
- Tests under `tests/` (parameterization suite + LSCM-system unit tests).

## Acceptance Criteria
- [ ] `RescaleMode { Off, LowStretchReference, GlobalAreaMin }` defined in
  `detail::lscm` and re-exported by both `AngleBasedLSCM` and `HierarchicalLSCM`.
- [ ] `set_rescale_mode(RescaleMode)` instance method + `RescaleMode`-accepting
  static `Compute` overloads on both solvers; default `Off`.
- [ ] On a mesh with a deliberately stretched auto-selected boundary edge,
  `GlobalAreaMin` yields total UV area ≈ total 3D area within tolerance, and
  beats `Off` on area preservation; `LowStretchReference` also improves over
  `Off`.
- [ ] `Off` reproduces current output exactly (all existing tests pass
  unchanged for both solvers).
- [ ] Rescale math has direct unit tests in the `detail::lscm` suite (factor
  computed correctly from known areas for both modes).
- [ ] Instance form (`set_rescale_mode` + `compute()`) covered for both solvers.
- [ ] Doxygen documents each mode, the origin-scaling pin caveat, and the
  planned 3.0 default change to `GlobalAreaMin`.
- [ ] Amalgamated single header regenerated and builds.

## Dependencies
- A7 (#42) — shared `detail::lscm` PinMap assembly path — merged.
- A8 (#64) — shared LSCM system-building extraction — merged.
- Coordinates with #96 (3.0 deprecation cleanup) for the eventual default flip.

## Out of Scope
- Smarter pin *selection/placement* (covered by closed A5/#2).
- Non-uniform / per-region rescaling or full SLIM-style area correction (F5).
- Changing the conformal solve itself; this is strictly a post-solve scale.

---
_Generated by Conductor. Review and edit as needed._
