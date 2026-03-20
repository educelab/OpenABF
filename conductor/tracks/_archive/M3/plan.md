# M3 Implementation Plan

## Phase 1 — Fix
Remove the `detail::erase_if` function from `HalfEdgeMesh.hpp`.

## Phase 2 — Verify
- Search codebase for any remaining references to `detail::erase_if` (should be none)
- Run `ctest`
