# M5 — split_edge and detail::filter cleanup

## GitHub Issues
- https://github.com/educelab/OpenABF/issues/76 — `detail::filter` has remove_if semantics
- https://github.com/educelab/OpenABF/issues/77 — split_edge: misleading exception messages and redundant vertex update
- https://github.com/educelab/OpenABF/issues/78 — split_edge: undetected non-manifold pinch when splitting interior fan edge at boundary vertex

## Summary
Three related code-quality issues surfaced during a review of
`HalfEdgeMesh::split_edge`. They all touch `include/OpenABF/HalfEdgeMesh.hpp`
within ~250 lines of each other and are bundled here to avoid churning the
same code three times.

1. **#76** — `detail::filter` (`HalfEdgeMesh.hpp:57-64`) is implemented with
   `std::remove_if`. It discards elements matching the predicate and keeps the
   rest — the opposite of the conventional `filter` semantic. Rename it to
   `detail::remove_if` (matches STL convention; minimum-churn fix).
2. **#77** — `split_edge` (`:1337-1342`, `:1360-1365`) throws
   `"No incoming/outgoing edges"` / `"Too many incoming/outgoing edges"` when
   the filtered list of *boundary* in/out edges has the wrong cardinality.
   Update messages to reflect that. Also remove the redundant
   `newFwd->prev->pair->vertex = newStart;` at `:1411` — it is set again at
   `:1415` in the boundary branch and is a no-op in the non-boundary branch.
3. **#78** — `split_edge` silently produces a non-manifold pinched vertex if
   the edge being split lives in the *middle* of a triangle fan at a boundary
   vertex (rather than at one end of the fan). Add a precondition check that
   throws a clear `MeshException` instead.

## Acceptance Criteria
- [ ] `detail::filter` is renamed to `detail::remove_if`; all call sites
      updated; predicate semantics unchanged.
- [ ] `split_edge` exception messages mention "boundary" so the failure mode
      is identifiable from the message alone.
- [ ] Redundant assignment at line 1411 is removed.
- [ ] `split_edge` throws a `MeshException` when called on an edge whose face
      is not adjacent to the existing boundary at either endpoint, instead of
      silently producing a non-manifold mesh.
- [ ] New unit test exercises the pinch case and expects the exception.
- [ ] All existing tests pass unchanged.
- [ ] `single_include/OpenABF/OpenABF.hpp` regenerated.

## Dependencies
None — independent of the active feature tracks.

## Out of scope
- F3 multi-component pipeline (separate track, separate PR).
- Larger refactors of `split_edge` (e.g. splitting the boundary/non-boundary
  branches into helpers) — defer to a follow-up if desired.
