# B9 Implementation Plan

Provide a recoverable mapping from a torn/parameterized mesh back to the
original input topology. Two unrecorded identity shifts must be invertible:
insert_face corner-order reversal, and split_edge seam-vertex duplication.
See spec.md and issue #100.

## Phase 1: Investigation & decision
- [ ] 1.1 Reproduce both shifts: (a) a mis-wound input face whose as-built corner
          order differs from input; (b) a torn seam whose duplicate vertex has no
          recorded link to its original — confirm neither is recoverable today
- [ ] 1.2 Survey re-insertion/duplication sites: insert_face/insert_faces (reversal),
          split_edge/split_path (duplication), clone_face_ (extraction) — confirm
          how each affects the original→result chain
- [ ] 1.3 Choose the approach (vertex origin tracking, returned remaps, per-face
          corner permutation, or a combination) and record the decision in spec.md

## Phase 2: Tests (write first)
- [ ] 2.1 Duplicate → original vertex recovery across split_path (incl. multi-seam)
- [ ] 2.2 Reversed-face corner-order recovery for a mis-wound input face
- [ ] 2.3 Full round-trip: packed/merged atlas corner → F2 maps → torn-mesh corner
          → B9 maps → original input face, input corner position, input vertex
- [ ] 2.4 Identity/no-op case: an untorn, correctly-wound mesh round-trips to identity

## Phase 3: Implementation
- [ ] 3.1 Implement duplicate → original vertex mapping in split_edge/split_path
- [ ] 3.2 Implement corner-order recovery for insert_face (reversed flag/permutation)
- [ ] 3.3 Ensure the mappings survive extract_connected_components (and compose with
          F2's vertex_map/face_map and vertex_source/face_source)
- [ ] 3.4 Document behavior + recovery API on insert_face/insert_faces/split_*
- [ ] 3.5 Regenerate single header (and update install list if a header is added)

## Phase 4: Verify
- [ ] 4.1 Run `ctest` — all suites pass
- [ ] 4.2 Run clang-format on changed files
- [ ] 4.3 Confirm single-header build compiles and runs
