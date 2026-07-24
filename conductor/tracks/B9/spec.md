# B9 — No recoverable mapping from torn seam-duplicate vertices back to the consumer's pre-HEM mesh

## GitHub Issue
https://github.com/educelab/OpenABF/issues/100

## Use case
Downstream consumers do **not** use `HalfEdgeMesh` as their primary mesh type.
They keep their own mesh (the "raw" / pre-HEM mesh — their own vertex array and
face list) and *wrap* it in a `HalfEdgeMesh` only to flatten: build the HEM with
`insert_vertices`/`insert_faces`, tear it (`split_path`), extract components,
parameterize, and pack/merge. They then build a **per-wedge UV map** from the
packed result and expect to apply it back onto their **pre-HEM mesh**, keyed by
their own `(face, corner)` identities.

## Summary
That round-trip is `atlas corner → (F2 maps) → torn-HEM corner → pre-HEM mesh`.
Three of the four hops already round-trip cleanly:

- **Face index.** `insert_faces` inserts faces in order
  (`include/OpenABF/HalfEdgeMesh.hpp:1040-1049`) and `split_path` only tears
  edges (it never re-inserts faces), so pre-HEM face `i` == HEM face `i`, and
  F2's `face_map`/`face_source` carry it back.
- **Corner order.** `insert_face` may auto-reverse a mis-wound face to keep the
  mesh manifold (`include/OpenABF/HalfEdgeMesh.hpp:1706-1748`), permuting a
  face's stored corner order relative to the raw input. This is **not** a
  blocker: the documented per-wedge recipe resolves corners by *vertex
  identity* against the consumer's own face, never by raw traversal index
  (`include/OpenABF/ChartPacking.hpp:108-123`), so the reversal is **absorbed**,
  not inverted. (See "Out of scope".)
- **Non-seam vertices.** `insert_vertices` preserves order and tearing only
  *appends* new vertices, so a non-duplicated HEM vertex keeps its pre-HEM
  index; F2's `vertex_map`/`vertex_source` carry it back.

The one hop that breaks is **seam-vertex identity**. Tearing duplicates seam
vertices via `insert_vertex(oldStart->pos)`
(`include/OpenABF/HalfEdgeMesh.hpp:1445,1465`), appending a new index whose only
link to its origin is a copied position. `split_path`/`split_edge` return `void`
and record no `duplicate → original` vertex map. So a torn-HEM corner that lands
on a seam duplicate **cannot be expressed in the consumer's pre-HEM vertex
namespace** — and identity-based corner resolution against the consumer's own
face (the very mechanism that absorbs the winding reversal) fails for exactly
those corners.

## Why this matters
The per-wedge UV recipe in `MultiChartFlatten.cpp` works entirely *within* the
torn HEM's own namespace, which is self-consistent. But the consumer's goal is
to land UVs on their **pre-HEM** mesh. For every non-seam corner that already
works; for a seam corner it silently cannot, because the duplicate has no
recoverable pre-HEM vertex identity. The same gap blocks scattering per-vertex
attributes captured on the pre-HEM mesh and emitting output indexed by the
consumer's original vertices. Discovered during F2 (PR #99).

## Goal
Make every torn/extracted vertex — original or seam duplicate — name its
**pre-HEM (input) vertex**, so the documented identity-keyed per-wedge recipe
resolves corners against the consumer's own faces for *all* corners, seams
included. Composed with F2's `vertex_map`/`face_map` and
`vertex_source`/`face_source`, a consumer can take any corner of a packed/merged
atlas and name the pre-HEM face, corner, and vertex it came from.

## Out of scope
- Removing auto-rewinding or seam duplication (both are intentional).
- Recording the `insert_face` corner-order permutation / a reversed flag.
  Identity-keyed corner resolution (`ChartPacking.hpp:108-123`) absorbs the
  reversal, so it does **not** need to be inverted for this use case. The corner
  order is recovered implicitly once seam-duplicate vertices carry their pre-HEM
  identity (this AC), by locating each vertex within the consumer's own face.

## Acceptance Criteria
- [ ] `split_edge`/`split_path` expose a recoverable **duplicate → pre-HEM**
      vertex mapping (e.g. an `origin` index stored on every vertex, set at
      construction and copied to duplicates, that survives extraction; or a
      returned/accumulated remap).
- [ ] The mapping survives `extract_connected_components` (`clone_face_` copies
      vertices, so an `origin` field rides along) and composes with F2's
      `vertex_map`/`vertex_source`.
- [ ] A worked path demonstrates the full round-trip: packed/merged atlas corner
      → (F2 maps) → torn-HEM corner → (B9 mapping) → pre-HEM face, corner, and
      vertex index — resolving the corner by vertex identity against the
      consumer's own face, with seam-duplicate corners resolving correctly.
- [ ] Unit tests on a mesh with at least one torn seam assert the duplicate →
      pre-HEM vertex mapping recovers the original identity, and that a seam
      corner of the packed atlas lands on the correct pre-HEM `(face, corner)`.
      Include a mis-wound input face to confirm the winding reversal is absorbed
      (the round-trip still lands the right corner without a corner-order record).
- [ ] An untorn, correctly-wound mesh round-trips to the identity mapping.
- [ ] `split_*` documentation describes the behavior and points to the recovery
      API; the `MultiChartFlatten.cpp` reference comment is updated to key
      corners against the pre-HEM mesh using the new vertex origin.
- [ ] Single-header regenerated; multiheader install list updated if a new
      header is introduced.

## Candidate approaches (decide in Phase 1)
1. **Vertex origin tracking (primary).** Store an `origin` (pre-HEM vertex index)
   set on construction and copied to duplicates by `split_edge`, so every vertex
   — original or duplicate — names its pre-HEM vertex. Survives
   `clone_face_`/extraction via the vertex copy path. With this, corners are
   located by identity in the consumer's own face and the winding reversal needs
   no separate record.
2. **Returned remaps.** `split_edge`/`split_path` return/accumulate
   `duplicate → original` pairs. Lighter-weight but does not survive extraction
   without the caller threading it through, and does not give a uniform
   "every vertex names its input vertex" accessor.

## Dependencies
- Independent of F2 (PR #99), but motivated by it; the F2 maps
  (`vertex_map`/`face_map`, `vertex_source`/`face_source`) are the downstream
  half of the chain B9 completes back to the pre-HEM mesh.
