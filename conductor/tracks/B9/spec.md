# B9 — insert_face auto-rewinding silently changes face corner order vs input

## GitHub Issue
https://github.com/educelab/OpenABF/issues/100

## Summary
`HalfEdgeMesh::insert_face` auto-reverses a mis-wound face to keep the mesh's
winding globally consistent (`include/OpenABF/HalfEdgeMesh.hpp:1706-1748`). The
reversal is desirable for making almost-manifold input manifold, but it
**silently changes a face's stored corner/vertex traversal order relative to the
raw input face list**, and the mesh records neither that a reversal occurred nor
the input→as-built corner permutation.

## Why this matters
Downstream consumers build a mesh from a known face list and later read back
per-face-corner data — applying a per-wedge UV map, transferring per-corner
attributes, emitting OBJ `vt`/`f` keyed by input corner order. For any reversed
face, the half-edge traversal order no longer matches the caller's input order,
with no signal. The discrepancy is silent.

Discovered during F2 (multi-chart UV packing, PR #99). The robust workaround is
to key per-wedge data by **vertex identity**, never by raw traversal index —
correct but non-obvious, and the trap is currently unguarded.

## Impact
- Silent mismatch: input corner order ≠ as-built corner order for reversed faces.
- Affects any round-trip of per-corner data through a `HalfEdgeMesh`.
- No API to detect a reversed face or recover the original corner order.

## Out of scope
- Removing the auto-rewinding behavior itself (it is intentional and useful).

## Acceptance Criteria
- [ ] Decision recorded on the chosen approach (see below).
- [ ] If a detection/recovery API is added: it reports, for each face, whether it
      was reversed at insertion and/or maps a corner index back to input order,
      with unit tests on a mesh containing at least one mis-wound input face.
- [ ] `insert_face` / `insert_faces` documentation prominently describes the
      auto-rewinding behavior and its effect on corner order.
- [ ] A test constructs a mesh with a deliberately mis-wound face and asserts the
      documented/observable behavior (and the recovery API if added).
- [ ] Single-header regenerated; multiheader install list updated if a new header
      is introduced (none expected).

## Candidate approaches (decide in Phase 1)
1. Record per-face reversal (flag) and/or store the input→as-built corner
   permutation, exposed via an accessor.
2. Provide a query mapping a corner index back to input order.
3. Opt-in "strict" insertion mode that throws on mis-wound input instead of
   silently reversing (caller fixes winding explicitly).
4. Minimum viable: document the behavior loudly + ship the identity-keying
   guidance (already drafted in the F2 MultiChartFlatten reference comment).

## Dependencies
- None. Independent of F2 (PR #99), though motivated by it.
