# A8 Implementation Plan

## Phase 1: Audit and design shared interface
1. Read `AngleBasedLSCM::ComputeImpl` and `HierarchicalLSCM::solveLSCMLevel` side by side
2. Identify exact shared lines and the parameters/return values needed
3. Design `detail::lscm::buildSystem<T, MeshPtr>(mesh, pin0, pin1)` → `SystemParts<T>`
   - `SystemParts` holds: `SparseMatrix<T> A`, `SparseMatrix<T> b` (or dense), `unordered_map<size_t,size_t> freeIdxTable`

## Phase 2: Extract utility (TDD)
1. Write `HLSCMInternal.BuildSystem_KnownMesh` test (pyramid, verify A/b dimensions and pin rows)
2. Implement `detail::lscm::buildSystem` in a new header `include/OpenABF/detail/LscmSystem.hpp`
   (or inline in `AngleBasedLSCM.hpp` at the bottom of the `detail` section)
3. Test passes

## Phase 3: Migrate AngleBasedLSCM
1. Replace duplicated logic in `AngleBasedLSCM::ComputeImpl` with call to `buildSystem`
2. All existing parameterization tests pass — numerical results identical

## Phase 4: Migrate HierarchicalLSCM
1. Replace duplicated logic in `solveLSCMLevel` with call to `buildSystem`
2. All HLSCM tests pass — numerical results identical

## Phase 5: Update amalgamation and finalize
1. Add new header to `single_include.json` if extracted to separate file
2. Run amalgamation script
3. Run clang-format
