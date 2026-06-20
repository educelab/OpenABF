# B9 Implementation Plan

Potential-bug track: insert_face auto-rewinding silently changes face corner
order vs the raw input, unrecorded. See spec.md and issue #100.

## Phase 1: Investigation & decision
- [ ] 1.1 Reproduce: build a mesh with a deliberately mis-wound face; confirm
          traversal corner order differs from input and nothing records it
- [ ] 1.2 Survey internal call sites that re-insert faces and could rewind
          (insert_faces, clone_face_ used by extract_connected_components, any
          merge/split helpers) and note where the discrepancy can leak
- [ ] 1.3 Choose an approach (flag / permutation accessor / strict mode / docs)
          and record the decision in spec.md

## Phase 2: Tests (write first)
- [ ] 2.1 Test asserting the chosen observable behavior on a mis-wound input face
- [ ] 2.2 If a recovery API is added: test that it maps corners back to input order
- [ ] 2.3 If a strict mode is added: test that mis-wound input throws under it

## Phase 3: Implementation
- [ ] 3.1 Implement the chosen approach
- [ ] 3.2 Document auto-rewinding + corner-order effect on insert_face/insert_faces
- [ ] 3.3 Regenerate single header (and update install list if a header is added)

## Phase 4: Verify
- [ ] 4.1 Run `ctest` — all suites pass
- [ ] 4.2 Run clang-format on changed files
- [ ] 4.3 Confirm single-header build compiles and runs
