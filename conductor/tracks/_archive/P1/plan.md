# P1 Implementation Plan

## Phase 1 — Tests
1. Write tests that exercise each iteration method in a range-based for loop:
   - `mesh->vertices()`, `mesh->edges()`, `mesh->faces()`
   - `mesh->vertices_interior()`, `mesh->vertices_boundary()`
   - `v->wheel()`, `mesh->outgoing_edges(idx)`
   - `mesh->boundaries()`
   - `face->edges()` (new named accessor, resolving A2)
2. These should serve as regression tests that the interface contract is preserved

## Phase 2 — Design range types
1. Create a lightweight `IterableView<Iterator>` or use `std::span` for contiguous storage
2. For linked-list traversals (`wheel`, `outgoing_edges`, `boundaries`, face edges),
   create a minimal input-iterator pair that wraps the linked-list walk
3. Keep the existing `std::vector`-returning methods as deprecated or remove them,
   deciding based on impact to examples and user-facing API

## Phase 3 — Migrate internal accessors
1. Update `HalfEdgeMesh` methods to return range views
2. Update `ABF::Compute`, `ABFPlusPlus::Compute`, `AngleBasedLSCM::Compute`,
   `InitializeAnglesAndWeights`, and all free functions to use the new iteration API
3. Update all tests and examples

## Phase 4 — Verify
- Run `ctest`
- Profile solver on a larger mesh to confirm allocation reduction
