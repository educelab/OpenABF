# A8 — Extract shared LSCM system-building logic

## GitHub Issue
https://github.com/educelab/OpenABF/issues/64

## Problem
`HierarchicalLSCM::solveLSCMLevel` duplicates approximately 150 lines of
`AngleBasedLSCM::ComputeImpl`. The duplicated sections are:

1. Pin placement and boundary-vertex handling
2. Free-vertex index table construction (`freeIdxTable`)
3. LSCM sparse system assembly (building the `A` matrix and `b` vector from edge angles)
4. UV extraction from the solver solution vector `x`

Any bug fix in one implementation must be manually applied to the other. The two have already
drifted (pin-axis selection, index table type). Future algorithm changes (e.g. A7 multi-pin
constraints) must be applied twice.

## Proposed design
Extract the shared logic into a `detail::lscm` namespace utility:

```cpp
namespace detail::lscm {
    // Build the LSCM sparse system for a mesh with two pinned vertices.
    // Returns (A, b, freeIdxTable) in canonical form.
    template <typename T, class MeshPtr>
    auto buildSystem(const MeshPtr& mesh, size_t pin0, size_t pin1)
        -> SystemParts<T>;  // { SparseMatrix A, SparseVector b, unordered_map freeIdx }
}
```

Both `AngleBasedLSCM::ComputeImpl` and `HierarchicalLSCM::solveLSCMLevel` call
`detail::lscm::buildSystem(...)` and handle only their solver-specific logic (warm-start,
solver dispatch, UV prolongation).

## Acceptance criteria
- [ ] `detail::lscm::buildSystem` (or equivalent) is a single implementation used by both solvers
- [ ] `solveLSCMLevel` and `AngleBasedLSCM::ComputeImpl` each call the shared utility
- [ ] All existing parameterization tests pass (numerical results unchanged)
- [ ] No new public API surface — `detail` namespace only
- [ ] The extracted utility is covered by at least one direct test (can use T5 infrastructure)
