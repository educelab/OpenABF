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
  shared frame; it does not own UV map types. F2 instead **documents** the
  reversal-robust per-wedge recipe (see Design Decisions) so consumers can
  build their own map correctly.
- **Writing UV coordinates back onto the input mesh's vertices.** The
  original `ParameterizeConnectedComponents` design from F3 was abandoned
  for this reason — a single mesh vertex can be shared across multiple
  charts after tearing, and any "write back to source" strategy silently
  destroys data for the shared-vertex case.

## Design Decisions (resolved 2026-06-20)

1. **Geometry-only, in-place.** `PackCharts` mutates each chart's
   `Vertex::pos` (translate + optional scale). Packing touches no topology
   and no vertex/face indices, so any `ExtractedComponent` back-maps the
   caller holds remain valid afterward.

2. **Per-wedge iteration is feasible and must key on vertex identity.**
   A `Face` stores no vertex list — corners come from walking the half-edge
   chain (`HalfEdgeMesh.hpp:408-474`), and `insert_face_` may **reverse a
   mis-wound face at insertion time** (`HalfEdgeMesh.hpp:1706-1748`), which
   permutes corner *position* relative to the raw input list. The HEM does
   not record that permutation. Therefore the robust wedge key is:

   ```
   wedge = (face_map[chart_face.idx], vertex_map[corner.vertex.idx])
   ```

   i.e. **original face id + original vertex id**, NOT (face, corner 0/1/2).
   Corner *vertex identity* is fully recoverable via the back-maps; corner
   *position* is not. This recipe is documented in the header.

3. **Scaling: absolute by default, opt-in normalize.** The LSCM default pin
   anchors one edge to its 3D length (`LSCMSystem.hpp:83-90`), so each chart
   arrives at an approximately consistent physical texel density. The 3D
   surface area itself is unavailable to `PackCharts` (parameterization
   overwrites `pos` with `{u,v,0}`, `AngleBasedLSCM.hpp:255-266`).
   - Default: **translate-only** packing — absolute scale untouched.
   - `normalize=true`: apply **one global uniform scale** to fit `[0,1]²`.
     A single global factor preserves relative sizes and cross-chart
     density; only absolute units change. (Per-chart rescaling is never
     done — it would destroy cross-chart density.)

4. **Layout: shelf packing.** Sort charts by bbox height descending, place
   left-to-right, wrap to a new shelf when the row exceeds the target width.
   Default target width = `sqrt(Σ chart bbox area)` (≈ square atlas);
   overridable via `PackOptions::target_width`.

5. **API shape.**
   ```cpp
   template <class T, class MeshType>
   struct PackOptions {
       bool normalize = false;            // fit packed atlas into [0,1]^2
       std::optional<T> target_width{};   // overrides sqrt-area heuristic
       T padding = T(0);                  // per-chart gutter, absolute units
   };
   struct PackResult { Vec<T,2> min, max; };   // packed atlas extent

   PackResult PackCharts(std::vector<typename MeshType::Pointer>& charts,
                         PackOptions opts = {});
   ```

6. **Degenerate input.** Empty list → return empty extent (no-op). Zero-area
   / single-point charts are placed by their (flat) bbox. Null pointer or a
   chart with zero vertices → `throw std::invalid_argument`.
   `static_assert(MeshType::Dim >= 2)`.

## Acceptance Criteria
- [ ] `PackCharts<MeshType>` free function accepts a list of parameterized
      mesh `Pointer`s and translates (and, when `normalize`, uniformly
      scales) each chart's vertex positions into a shared coordinate frame.
- [ ] Default mode preserves absolute scale (translate-only); cross-chart
      relative sizes are never distorted.
- [ ] `normalize=true` fits the packed atlas into `[0,1]²` via a single
      global uniform scale.
- [ ] Shelf-packing with the ~square target-width heuristic, overridable.
- [ ] No charts' bounding boxes overlap (padding respected); the packed set
      is contained in the returned extent (and in `[0,1]²` when normalized).
- [ ] Edge cases handled per Design Decision 6.
- [ ] Header documents the vertex-identity per-wedge recipe (Decision 2).
- [ ] Tests: synthetic 2D charts for deterministic geometric assertions plus
      one end-to-end `tear → extract → LSCM → pack` test that exercises
      wedge recovery via the back-maps.
- [ ] Documented with complexity notes.

## Dependencies
- F3 (multi-component extraction via
  `HEM::extract_connected_components()`) is **complete** (PRs #80, #81) and
  provides `ExtractedComponent{mesh, vertex_map, face_map}` — the natural
  input to this function and the source of the back-maps used by the
  documented wedge recipe.
