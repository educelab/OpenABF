# A7 Implementation Plan

Branch: `a7-multi-pin-lscm`.

Scope expanded 2026-06-14: multi-pin now covers both `AngleBasedLSCM` and
`HierarchicalLSCM`. Most of the work lives in a shared `detail::lscm`
assembly helper used by both solvers.

## Phase 1: Tests (Red) — done
- [x] 1.1 ABLSCM: 3-pin pyramid test — pin three vertices with explicit UV
  coordinates; verify each lands exactly at its specified UV.
- [x] 1.2 ABLSCM: instance-form test — `setPins(PinMap)` + `compute()` produces
  the same result as the static `Compute(mesh, PinMap)` overload.
- [x] 1.3 HLSCM: 3-pin pyramid test (single-level fallback path) — verify each
  pin lands at its specified UV.
- [x] 1.4 HLSCM: 3-pin hemisphere test (multi-level path) — verify pins survive
  every hierarchy level and land at the specified UVs.
- [x] 1.5 HLSCM: instance-form test — `setPins(PinMap)` + `compute()`.

Confirmed Red: `cmake --build . --target OpenABF_TestParameterization` fails
with 12 unknown-symbol errors (`PinMap`, `setPins`); no regressions in
existing tests.

## Phase 2: Migrate `detail::lscm::buildSystem` to PinMap-only
Per user direction 2026-06-14: no two-arg compat. PinMap is the only explicit
pin API in this codebase after this migration.

- [x] 2.1 Redesign `buildSystem` to a single signature
  `buildSystem(mesh, PinMap)` (no axis-snap inside; caller provides UVs
  verbatim). Mutates each pinned vertex's `pos` to its specified UV, then
  emits the LSCM least-squares system with those values in `bFixed`.
- [x] 2.2 Migrate `LSCMSystemBuild.*` tests + `BuildPyramidSystem` helper to
  the new signature.

## Phase 3: ABLSCM full migration to PinMap (with deprecated 2-pin shim)
Per user direction 2026-06-14: retain `Compute(mesh, p0, p1)` and
`setPinnedVertices(p0, p1)` as `[[deprecated]]` shims slated for removal in
3.0. Both forward to the PinMap path with axis-snap UVs.

- [x] 3.1 Define `using PinMap = std::vector<std::pair<std::size_t, Vec<T,2>>>`
  in `AngleBasedLSCM`.
- [x] 3.2 Add private helpers `selectBoundaryPair(mesh)` and
  `autoPlacePair(mesh, p0Idx, p1Idx) -> PinMap` (the existing axis-snap
  convention). `autoSelectPins(mesh) = autoPlacePair(mesh,
  selectBoundaryPair(mesh)...)`.
- [x] 3.3 Rewrite `Compute(mesh)`: build PinMap via `autoSelectPins`, then
  call `ComputeImpl(mesh, pins)`.
- [x] 3.4 Add `Compute(mesh, PinMap)` static overload (validates pins ≥ 2,
  unique indices, in-range).
- [x] 3.5 Add `setPins(PinMap)` and `pins_` (`std::optional<PinMap>`).
- [x] 3.6 Restore `Compute(mesh, p0Idx, p1Idx)` and
  `setPinnedVertices(p0, p1)` as `[[deprecated]]` shims that build an
  axis-snap PinMap and delegate. Update `compute()` dispatch to handle
  legacy indices stored by the shim.
- [x] 3.7 Keep `AngleBasedLSCM_ExplicitPins`, `_Reversed`,
  `_SetPinnedVertices` — they now exercise the deprecated-shim path and
  must continue to pass until 3.0.

## Phase 4: HLSCM full migration to PinMap (with deprecated 2-pin shim)
- [x] 4.1 `DecimationMesh::build` takes `std::vector<std::size_t>` pinIndices;
  flag each as `isPinned_[idx] = true`.
- [x] 4.2 `buildHierarchy` takes the pinIndices vector and forwards.
- [x] 4.3 `solveLSCMLevel` takes a PinMap; map each pin's original idx to its
  level-local idx, write the user's UV into the level mesh, run buildSystem,
  copy each pin back into the output UV vector.
- [x] 4.4 Update initial-guess builder to skip every pin (not just p0/p1).
- [x] 4.5 In `HierarchicalLSCM`: add `PinMap` alias, `autoSelectPins` helper
  using the existing axis-snap convention, `Compute(mesh)`,
  `Compute(mesh, PinMap)`, `setPins`, `pins_`, and rewrite `compute()`
  dispatch.
- [x] 4.6 Single-level fallback delegates to
  `AngleBasedLSCM::Compute(mesh, PinMap)`.
- [x] 4.7 Restore `Compute(mesh, pin0Idx, pin1Idx)` and
  `setPinnedVertices(p0, p1)` as `[[deprecated]]` shims (parallel to ABLSCM).
- [x] 4.8 Migrate internal tests that pass the two-arg form to the new
  signature (`std::vector<std::size_t>` for `DecimationMesh::build` and
  `buildHierarchy`; `PinMap` for `solveLSCMLevel`).
- [x] 4.9 Keep `HLSCM.ExplicitPins`, `HLSCM.SetPinnedVertices` — they now
  exercise the deprecated-shim path and must continue to pass until 3.0.

## Phase 5: Verify — done
- [x] 5.1 Run `ctest` — 6/6 binaries, 50/50 parameterization assertions pass.
- [x] 5.2 New multi-pin tests pass (5/5 new + 5 deprecated-shim tests still
  green).
- [x] 5.3 Run `git clang-format`; regenerate amalgamated header.

## Phase 6: Conductor & GitHub
- [x] 6.1 PR #93 opened against #42; PR description notes expanded scope
  (HLSCM included) and the deprecated-shim plan for 3.0 removal.
- [x] 6.2 Update `tracks.md` (move A7 to archived once PR lands).
