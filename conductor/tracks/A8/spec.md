# A8 — Extract shared LSCM system-building logic

## GitHub Issue
https://github.com/educelab/OpenABF/issues/64

## Downstream
- [P6](../P6/spec.md) (#66) is blocked on A8 — its sin/cos precomputation
  targets the face-assembly loop being extracted here.
- P4 and P5 are mostly orthogonal and can run in parallel after A8 lands.

## Problem
`HierarchicalLSCM::solveLSCMLevel` duplicates approximately 150 lines of
`AngleBasedLSCM::ComputeImpl`. The duplicated sections are:

1. Pin placement on UV axes
2. Free-vertex index table construction (`freeIdxTable`)
3. LSCM sparse system assembly (building the `A` matrix and `b` vector from edge angles)
4. UV extraction from the solver solution vector `x`

Any bug fix in one implementation must be manually applied to the other. The two
have already drifted: `freeIdxTable` is `std::map` in `AngleBasedLSCM` and
`std::unordered_map` in `HierarchicalLSCM`. Future algorithm changes
(e.g. A7 multi-pin constraints) must be applied twice.

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
- [ ] All existing parameterization tests pass (numerical results unchanged, bit-identical on reference meshes)
- [ ] No new public API surface — `detail` namespace only
- [ ] `freeIdxTable` container type is unified to `std::unordered_map<std::size_t, std::size_t>` (O(1) lookup; ABLSCM only uses `.at()` so the change is numerically inert)
- [ ] The extracted utility is covered by direct tests for (a) dimensions of `A`/`b`, (b) free-vertex index table population, (c) pin-row contributions land in `b` not `A`
