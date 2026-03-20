# A5 Implementation Plan

## Phase 1 — Tests
1. Add a test that runs LSCM with a manually specified pin pair and verifies the
   pinned vertices end up at the expected positions
2. Add a test using the longest-boundary-edge heuristic

## Phase 2 — API change
1. Add an overload of `Compute` that accepts explicit pinned vertex indices:
   ```cpp
   static void Compute(typename Mesh::Pointer& mesh,
                       std::size_t pin0, std::size_t pin1);
   ```
2. Extract the pinned-vertex selection logic from `Compute` into a separate
   `SelectPinnedEdge` free function (longest-boundary-edge variant as optional)
3. Refactor the existing `Compute` to use `SelectPinnedEdge` with the default heuristic
4. Add a `setPinnedVertices(std::size_t, std::size_t)` method on the class form

## Phase 3 — Verify
- Run `ctest`
- Check that the default pyramid test still produces the same result
