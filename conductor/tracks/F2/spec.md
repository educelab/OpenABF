# F2 — Multi-chart UV packing

## GitHub Issue
https://github.com/educelab/OpenABF/issues/18

## Summary
After splitting a mesh and parameterizing each connected component, there is no
utility to pack the resulting UV charts into a [0,1]² atlas. This is a stated
project goal.

## Acceptance Criteria
- [ ] A `PackCharts` (or similar) free function accepts a list of parameterized
  `HalfEdgeMesh` instances and packs their UV coordinates into a shared atlas space
- [ ] Charts are scaled uniformly (preserving relative area ratios) before packing
- [ ] At minimum, a shelf-packing or simple rectangle-packing algorithm is implemented
- [ ] Output is either UV coordinate updates in-place or a separate UV map structure
- [ ] Tests verify no charts overlap and all fit within [0,1]²
- [ ] Documented with complexity notes

## Dependencies
- F3 (multi-component extraction) provides the natural input to this function;
  implement F3 first so the end-to-end pipeline can be tested together
