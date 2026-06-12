# F2 — Multi-chart packing in a common world coordinate frame

## GitHub Issue
https://github.com/educelab/OpenABF/issues/18

## Summary
After tearing a mesh with `split_path` and parameterizing each connected
component independently (see F3), each chart is a free-floating 2D mesh in
its own coordinate frame. Packing transforms those flattened sub-meshes into
a **shared world coordinate frame** — typically scaled to preserve relative
area ratios and laid out so they do not overlap.

## Scope
F2 operates on geometry, not on UV map structures. The input is a list of
parameterized `HalfEdgeMesh` instances; the output rewrites each sub-mesh's
2D vertex positions so the whole set lives in a common frame. Downstream
consumers (mesh IO, atlas writers, packers) can read the resulting positions
and emit per-corner `vt` entries, build atlases, etc., from there.

## Explicitly out of scope
- **A `UVMap` / per-wedge UV map class.** UV-map data structures and
  per-wedge bookkeeping belong in a separate repository (consumer-side).
  OpenABF's job is to produce parameterized meshes and place them in a
  shared frame; it does not own UV map types.
- **Writing UV coordinates back onto the input mesh's vertices.** The
  original `ParameterizeConnectedComponents` design from F3 was abandoned
  for this reason — a single mesh vertex can be shared across multiple
  charts after tearing, and any "write back to source" strategy silently
  destroys data for the shared-vertex case.

## Acceptance Criteria
- [ ] `PackCharts` (or similar) free function accepts a list of
      parameterized `HalfEdgeMesh` instances and translates / scales each
      chart's vertex positions so the whole set fits into a common world
      coordinate frame.
- [ ] Charts are scaled uniformly across the set to preserve relative area
      ratios.
- [ ] At minimum, a shelf-packing or simple rectangle-packing algorithm is
      implemented.
- [ ] No charts overlap; all charts fit within the target frame.
- [ ] Tests verify non-overlap and target-frame containment.
- [ ] Documented with complexity notes.

## Dependencies
- F3 (multi-component extraction via
  `HEM::extract_connected_components()`) provides the natural input to this
  function. Implement F3 first so the end-to-end pipeline can be tested
  together.
