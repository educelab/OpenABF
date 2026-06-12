# M5 Implementation Plan

## Phase 1 — Tests
1. Add a unit test in `tests/src/TestHalfEdgeMesh.cpp` that constructs a
   boundary vertex with a 3-triangle fan and attempts to split an edge in the
   middle of the fan (the pinch case from #78). Expect `MeshException`.
2. Confirm the test FAILS against current `split_edge` (silent corruption
   rather than an exception).

## Phase 2 — Implementation
1. **#76** — Rename `detail::filter` → `detail::remove_if` at
   `include/OpenABF/HalfEdgeMesh.hpp:59`. Update all 6 call sites
   (`:1041, :1054, :1333, :1335, :1356, :1358`). Predicates unchanged.
2. **#77** — Replace the two exception messages in `split_edge` (at the new
   line numbers after rename) with `"Non-manifold boundary vertex at split
   endpoint"`. Drop the `== 0` branch entirely since it is unreachable given
   the surrounding `startOnBoundary`/`endOnBoundary` checks. Remove
   `newFwd->prev->pair->vertex = newStart;` at line 1411.
3. **#78** — Before mutating in `split_edge`, validate the precondition: when
   `startOnBoundary`, walk from `startOut->pair` via `prev` (or from `startIn`
   via `next`) and confirm the traversal reaches `oldFwd` without leaving
   `oldFwd`'s face. If not, throw `MeshException("split_edge: would create
   non-manifold pinched vertex")`. Symmetric check at `oldEnd`.

## Phase 3 — Verify
- `git clang-format` against staged changes.
- Run `python3 thirdparty/amalgamate/amalgamate.py -c single_include.json -s .`
  to regenerate `single_include/OpenABF/OpenABF.hpp`.
- Run full `ctest` and confirm all tests pass.

## Phase 4 — Ship
- Push branch `cleanup/split-edge-and-filter`.
- Open draft PR with `Fixes #76`, `Fixes #77`, `Fixes #78` in the body so the
  Development feature auto-links all three issues.
