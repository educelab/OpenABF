# A7 Implementation Plan

## Phase 1: Tests
- [ ] 1.1 Add a test pinning three vertices of the pyramid with explicit UV coordinates and verifying each lands at the specified position
- [ ] 1.2 Add a test using the instance form `setPins()` + `compute()`

## Phase 2: Implementation
- [ ] 2.1 Define `using PinMap = std::vector<std::pair<std::size_t, Vec<T, 2>>>` inside `AngleBasedLSCM`
- [ ] 2.2 Refactor `ComputeImpl` (introduced by A5) to accept a `PinMap` instead of two vertex pointers; build `bFixed` and the free/fixed split from the map
- [ ] 2.3 Add `Compute(mesh, PinMap)` static overload
- [ ] 2.4 Add `setPins(PinMap)` instance method and `pins_` member (`std::optional<PinMap>`)
- [ ] 2.5 Update `compute()` to call the PinMap overload when `pins_` is set

## Phase 3: Verify
- [ ] 3.1 Run `ctest`
- [ ] 3.2 Confirm all existing tests (default and two-pin paths) pass unchanged
