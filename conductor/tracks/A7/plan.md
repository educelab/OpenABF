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

## Phase 2: Shared assembly helper
- [ ] 2.1 Add `detail::lscm::buildSystemMultiPin(mesh, pins)` in
  `LSCMSystem.hpp`. Takes a `std::vector<std::pair<size_t, Vec<T,2>>>`. No
  axis-snap; pin positions are written verbatim to the mesh and into `bFixed`.
  Returns the existing `SystemParts<T>` shape.
- [ ] 2.2 Keep the two-pin `buildSystem(mesh, p0, p1)` overload — it stays the
  home of the LSCM axis-snap auto-placement convention used by the
  default/two-pin code paths.

## Phase 3: ABLSCM multi-pin
- [ ] 3.1 Define `using PinMap = std::vector<std::pair<std::size_t, Vec<T,2>>>`
  inside `AngleBasedLSCM`.
- [ ] 3.2 Add `Compute(mesh, PinMap)` static overload. Validates pins ≥ 2,
  unique indices, in-range; dispatches through `buildSystemMultiPin`.
- [ ] 3.3 Add `setPins(PinMap)` and `pins_` (`std::optional<PinMap>`).
- [ ] 3.4 Update `compute()` to prefer `pins_` over `pinnedVertices_` when set.

## Phase 4: HLSCM multi-pin
- [ ] 4.1 Generalize `DecimationMesh::build` to accept a span/vector of pin
  indices and flag each as `isPinned_[i] = true`.
- [ ] 4.2 Generalize `buildHierarchy` to forward the pin set to `build`.
- [ ] 4.3 Generalize `solveLSCMLevel` to take a PinMap; map each pin's original
  idx to its level-local idx, write the user's UV into the level-mesh vertex
  positions, build via `buildSystemMultiPin`, and write each pin back in the
  output UV vector.
- [ ] 4.4 Update initial-guess builder to skip *every* pin (not just p0/p1).
- [ ] 4.5 Add `HierarchicalLSCM::Compute(mesh, PinMap)` and `setPins(PinMap)`.
- [ ] 4.6 Update `compute()` to prefer `pins_` over `pinnedVertices_`.

## Phase 5: Verify
- [ ] 5.1 Run `ctest` — all existing tests pass unchanged.
- [ ] 5.2 New multi-pin tests pass.
- [ ] 5.3 Run `git clang-format`; regenerate amalgamated header.

## Phase 6: Conductor & GitHub
- [ ] 6.1 Open PR against issue #42; note expanded scope (HLSCM included) in PR
  description.
- [ ] 6.2 Update `tracks.md` (move A7 to archived once PR lands).
