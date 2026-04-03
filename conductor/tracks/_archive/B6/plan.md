# B6 Implementation Plan

## Phase 1 — Test
Add a test in `TestMeshIO.cpp` that attempts to read a PLY string/file with a
missing `z` property and asserts a `std::runtime_error` is thrown.

## Phase 2 — Fix
After the `vmap` population loop, add validation:
```cpp
std::array<bool, 3> vmapFound{false, false, false};
// set to true when x, y, z are found respectively
// after loop:
if (!vmapFound[0] || !vmapFound[1] || !vmapFound[2]) {
    throw std::runtime_error("PLY vertex element missing required x/y/z properties");
}
```

## Phase 3 — Verify
- Run `ctest`
