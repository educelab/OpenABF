# B9 — No recoverable mapping from a torn/parameterized mesh back to the original input topology

## GitHub Issue
https://github.com/educelab/OpenABF/issues/100

## Summary
A consumer who builds a `HalfEdgeMesh` from a known face list, tears it
(`split_path`), extracts components, parameterizes, and packs them has **no
recorded path back to their original input topology**. Two distinct,
*independently unrecorded* identity shifts sit between the input and the result:

1. **insert_face winding reversal (corner order).** `insert_face` auto-reverses
   a mis-wound face to keep the mesh globally manifold
   (`include/OpenABF/HalfEdgeMesh.hpp:1706-1748`). This is intended and useful,
   but it permutes a face's stored corner order relative to the raw input face
   list, and the permutation is not recorded.
2. **split_edge vertex duplication (vertex identity).** Tearing duplicates seam
   vertices via `insert_vertex(oldStart->pos)`
   (`include/OpenABF/HalfEdgeMesh.hpp:1445,1465`), appending new indices.
   `split_path`/`split_edge` return `void` and record no `duplicate → original`
   vertex map; the only thing the duplicate carries is a copied position.

These **compose**: rewinding happens at build time, then seam splitting changes
which vertex identity sits at a (still position-stable) corner. Working purely
in the post-split mesh's own namespace is fine and single-level (this is what
the F2 `PackCharts`/`MultiChartFlatten` per-wedge recipe does — key by position,
resolve identity within M′). But the moment a user needs to relate results back
to **their original input** — original vertex indices and original (input) face
corner order — both shifts must be inverted, and neither is currently
recoverable.

## Why this matters
Round-tripping per-corner / per-vertex data to the *caller's own topology* is a
normal need: applying an externally authored per-wedge UV map, transferring
per-vertex attributes captured before tearing, or emitting output indexed the
way the caller supplied geometry. Today that bridge is silently impossible for
any reversed face or any seam-duplicated vertex. Discovered during F2 (PR #99).

## Goal
Provide a **recoverable mapping from the torn/extracted result back to the
original input mesh**, covering both shifts:
- duplicate seam vertex → original input vertex, and
- as-built face corner order → raw input face corner order.

Composed with the existing `ExtractedComponent::vertex_map`/`face_map` (and
`MergedMesh::vertex_source`/`face_source` from F2), a consumer should be able to
take any corner of a packed/merged atlas and name the original input mesh's
face, input corner position, and input vertex it came from.

## Out of scope
- Removing auto-rewinding or seam duplication (both are intentional).

## Acceptance Criteria
- [ ] `split_edge`/`split_path` expose a recoverable **duplicate → original**
      vertex mapping (e.g. returned map, accumulator, or an origin index stored
      on duplicated vertices that survives extraction).
- [ ] `insert_face`/`insert_faces` expose whether a face was reversed and/or a
      way to recover the raw input corner order (reversed flag or
      input→as-built corner permutation accessor).
- [ ] A worked path demonstrates full round-trip: packed/merged atlas corner →
      (F2 maps) → torn-mesh corner → (B9 maps) → original input face, input
      corner position, and input vertex index.
- [ ] Unit tests on a mesh with (a) at least one deliberately mis-wound input
      face and (b) at least one torn seam, asserting both inverse mappings
      recover the original identities.
- [ ] `insert_face`/`split_*` documentation describes the behavior and points to
      the recovery API.
- [ ] Single-header regenerated; multiheader install list updated if a new
      header is introduced.

## Candidate approaches (decide in Phase 1)
1. **Vertex origin tracking.** Store an `origin` (original input vertex index)
   that is set on construction and copied to duplicates by `split_edge`, so
   every vertex — original or duplicate — names its input vertex. Survives
   `clone_face_`/extraction via the vertex copy path.
2. **Returned remaps.** `split_edge`/`split_path` return/accumulate
   `duplicate → original` pairs; a separate per-face reversal record handles
   corner order.
3. **Per-face corner permutation.** Record, per face, the input→as-built corner
   order (a reversed flag suffices for triangles; a rotation+reversal for
   general polygons), with an accessor.
4. **Combination + docs.** Likely (1)+(3): identity via vertex origin, corner
   order via per-face reversal record, plus prominent documentation and the
   identity-keying guidance already drafted in the F2 reference comment.

## Dependencies
- Independent of F2 (PR #99), but motivated by it; the F2 maps
  (`vertex_map`/`face_map`, `vertex_source`/`face_source`) are the downstream
  half of the chain B9 completes back to the input.
