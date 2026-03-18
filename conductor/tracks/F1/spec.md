# F1 — HLSCM (Hierarchical LSCM)

## GitHub Issue
https://github.com/educelab/OpenABF/issues/5

## Summary
No hierarchical LSCM implementation exists. The approach described in the reference
paper provides better parameterizations in high-curvature regions by hierarchically
building and refining the parameterization.

## Related
GitHub issue #5
Reference: https://static.aminer.org/pdf/PDF/000/593/591/least_squares_conformal_maps_for_automatic_texture_atlas_generation.pdf

## Acceptance Criteria
- [x] `HLSCM<T>` class is implemented with `compute()` and static `Compute()` API
  matching the style of `AngleBasedLSCM`
- [x] Works with the same `HalfEdgeMesh` types as `AngleBasedLSCM`
- [x] Tests verify correctness against known mesh configurations
- [x] Documented with algorithm citation
- [x] Included in `OpenABF.hpp`

## Dependencies
- A5 (pinned edge selection) may inform HLSCM's pinning strategy
