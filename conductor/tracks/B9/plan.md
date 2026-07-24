# B9 Implementation Plan

Make every torn/extracted vertex name its pre-HEM (input) vertex, so a consumer
who wraps their own mesh in a HalfEdgeMesh to flatten can build a per-wedge UV
map that lands back on their pre-HEM mesh — seam-duplicated corners included.
The `insert_face` winding reversal is absorbed by identity-keyed corner
resolution and needs no separate record (see spec.md and issue #100).

## Phase 1: Investigation & decision
- [ ] 1.1 Reproduce the gap: tear a seam, confirm the duplicate vertex has no
          recorded link to its pre-HEM origin and that an atlas seam corner
          cannot be resolved to the consumer's `(face, corner)` today
- [ ] 1.2 Confirm the already-recoverable hops to keep scope tight: face index
          (insert_faces order + split_path never re-inserts), non-seam vertex
          index (insert_vertices order + tearing only appends), and that
          identity-keyed corner resolution absorbs the winding reversal
- [ ] 1.3 Survey duplication/copy sites: split_edge/split_path (duplication) and
          clone_face_ (extraction vertex copy) — confirm an `origin` field rides
          the copy path; choose vertex-origin tracking vs returned remaps and
          record the decision in spec.md

## Phase 2: Tests (write first)
- [ ] 2.1 Duplicate → pre-HEM vertex recovery across split_path (incl. multi-seam)
- [ ] 2.2 Full round-trip: packed/merged atlas corner → F2 maps → torn-HEM corner
          → B9 mapping → pre-HEM face, corner, vertex; assert a seam corner lands
          on the correct pre-HEM (face, corner) by identity
- [ ] 2.3 Mis-wound input face: confirm the winding reversal is absorbed — the
          round-trip lands the correct corner with no corner-order record
- [ ] 2.4 Identity/no-op case: an untorn, correctly-wound mesh round-trips to the
          identity mapping

## Phase 3: Implementation
- [ ] 3.1 Implement duplicate → pre-HEM vertex mapping in split_edge/split_path
          (origin set at construction, copied to duplicates)
- [ ] 3.2 Ensure the mapping survives extract_connected_components (clone_face_
          vertex copy) and composes with F2's vertex_map/vertex_source
- [ ] 3.3 Document behavior + recovery API on split_*; update the
          MultiChartFlatten.cpp reference comment to key corners against the
          pre-HEM mesh using the vertex origin
- [ ] 3.4 Regenerate single header (and update install list if a header is added)

## Phase 4: Verify
- [ ] 4.1 Run `ctest` — all suites pass
- [ ] 4.2 Run clang-format on changed files
- [ ] 4.3 Confirm single-header build compiles and runs
