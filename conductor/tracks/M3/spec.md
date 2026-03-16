# M3 — `detail::erase_if` name misleads callers

## GitHub Issue
https://github.com/educelab/OpenABF/issues/27

## Summary
`detail::erase_if` takes its container argument by value, modifies the copy, and
returns it — making it a *filtered copy* utility. This conflicts with the name
`erase_if`, which implies in-place mutation matching C++20's `std::erase_if`. The
function has 6 active call sites in `HalfEdgeMesh.hpp` (lines 800, 813, 1094,
1096, 1117, 1119), all of which correctly capture the return value.

## Accepted Approach
Rename `detail::erase_if` → `detail::filter`. The library targets C++17, so
`std::erase_if` (C++20) and `std::ranges::filter_view` are unavailable; keeping
a named helper is appropriate. The rename makes the copy semantics self-evident.

## Acceptance Criteria
- [ ] `detail::erase_if` renamed to `detail::filter` in `HalfEdgeMesh.hpp`
- [ ] All 6 call sites updated
- [ ] All tests pass
- [ ] No change in behavior at any call site

## Dependencies
None.
