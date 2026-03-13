# A6 Implementation Plan

## Phase 1 — Refactor
1. Extract a helper lambda inside `Compute`:
   ```cpp
   auto addContrib = [&](std::size_t row, const EdgePtr& e, T c, T s) {
       if (e->vertex == p0) {
           // push to tripletsB (fixed)
       } else if (e->vertex == p1) {
           // push to tripletsB (fixed)
       } else {
           auto freeIdx = freeIdxTable.at(e->vertex->idx);
           // push to tripletsA (free)
       }
   };
   ```
2. Replace the three repeated blocks (for e0, e1, e2) with calls to `addContrib`
3. Ensure the mathematical layout matches equations in Lévy et al. (2002) in comments

## Phase 2 — Verify
- Run `ctest`
- All parameterization tests must produce bit-identical results
