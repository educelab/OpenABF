# F1 Implementation Plan

## Phase 1 — Research & Spec
1. Read the Lévy et al. paper section on hierarchical LSCM
2. Document the algorithm steps in a design note inside this track directory
3. Identify any additional data structures or mesh traits needed

## Phase 2 — Tests
1. Write failing tests for `HLSCM::Compute` on the existing pyramid fixture
2. Write tests for a more complex mesh where HLSCM is expected to outperform LSCM

## Phase 3 — Implementation
1. Create `include/OpenABF/HLSCM.hpp`
2. Implement `HLSCM<T, MeshType, Solver>` following the `AngleBasedLSCM` pattern
3. Add to `include/OpenABF/OpenABF.hpp`
4. Update `single_include.json` for the amalgamated header

## Phase 4 — Verify
- Run `ctest`
- Compare HLSCM vs LSCM output on complex meshes
- Build docs and review Doxygen output
