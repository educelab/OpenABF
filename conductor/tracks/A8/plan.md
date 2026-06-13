# A8 Implementation Plan

## Phase 1: Audit and design shared interface
- [x] 1.1 Read `AngleBasedLSCM::ComputeImpl` and `HierarchicalLSCM::solveLSCMLevel` side by side; confirm pin-axis placement and `addContrib` lambda are already in sync — the only remaining drift is the `freeIdxTable` container type
- [x] 1.2 Record the design decisions below in the plan (these are deliberate up-front choices, not Phase-3 discoveries):
  - **Signature**: `detail::lscm::buildSystem<T, MeshType>(meshPtr, p0VertPtr, p1VertPtr) -> SystemParts<T>` where `SystemParts<T>` holds `SparseMatrix<T> A`, `SparseMatrix<T> b`, `std::unordered_map<std::size_t, std::size_t> freeIdxTable`
  - **Responsibilities**: buildSystem performs pin placement on UV axes (mutates `p0->pos`, `p1->pos`) AND assembly. Pin *selection* (boundary-search or index→VertPtr) stays in the callers — that is where the two solvers correctly differ
  - **`freeIdxTable` type**: `std::unordered_map<std::size_t, std::size_t>` in the unified utility. ABLSCM only uses `.at()` lookup, so the swap from `std::map` is numerically inert
  - **`addContrib` lambda**: stays as a local lambda inside `buildSystem` — not exposed; revisit at A7 if multi-pin needs composability
  - **Header location**: new `include/OpenABF/detail/LSCMSystem.hpp`. Both `AngleBasedLSCM.hpp` and `HierarchicalLSCM.hpp` include it. HLSCM keeps its existing include of `AngleBasedLSCM.hpp` for `detail::SolveLeastSquares` / `detail::is_instance_of_v` (out of scope for A8)

## Phase 2: Extract utility (TDD)
- [x] 2.1a Write `LSCMSystemBuild.Dimensions_KnownMesh` test (pyramid) — assert `A` is `2·numFaces × 2·numFree`, `b` is `2·numFaces × 1`
- [x] 2.1b Write `LSCMSystemBuild.FreeIdxTable_Population` test — assert size = `numVerts - 2`, contains all non-pin vertex indices, no pin indices
- [x] 2.1c Write `LSCMSystemBuild.PinRowsLandInB` test — assert pin-row contributions appear in `b` (via the bFixed contraction) and not in `A`
- [x] 2.2 Implement `detail::lscm::buildSystem` in new header `include/OpenABF/detail/LSCMSystem.hpp` (transitively reached via `AngleBasedLSCM.hpp`; no `single_include.json` edit needed — amalgamator follows the include graph)
- [x] 2.3 All three new tests pass; full ctest suite (6/6) still passes

## Phase 3: Migrate AngleBasedLSCM
- [x] 3.1 Replace duplicated logic in `AngleBasedLSCM::ComputeImpl` with call to `detail::lscm::buildSystem`
- [x] 3.2 All existing parameterization tests pass — `EXPECT_FLOAT_EQ` (`Parameterization.AngledBasedLSCM`) and `EXPECT_DOUBLE_EQ` (`Parameterizations.AngleBasedLSCM_Double`) against hardcoded baseline UVs both pass → bit-identical within 4 ULP

## Phase 4: Migrate HierarchicalLSCM
- [x] 4.1 Replace duplicated logic in `solveLSCMLevel` with call to `detail::lscm::buildSystem`
- [x] 4.2 All HLSCM tests pass — including `HLSCM.MultiLevelHierarchy`, `HLSCM.ABFPlusPlusAnglePreservation` (T7 curved mesh), and `HLSCMInternal.SolveLSCMLevel_KnownMesh` which directly verifies parity with `AngleBasedLSCM`

## Phase 5: Amalgamation and finalize
- [ ] 5.1 Verify amalgamation script picks up new header (`single_include.json` updated in 2.2); run `python3 thirdparty/amalgamate/amalgamate.py -c single_include.json -s .`
- [ ] 5.2 Run `git clang-format` and re-stage
