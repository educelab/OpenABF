# M3 — `detail::erase_if` name misleads callers

## GitHub Issue
https://github.com/educelab/OpenABF/issues/27

## Summary
`detail::erase_if` takes its container argument by value, modifies the copy, and
returns it — making it a *filtered copy* utility. This conflicts with the name
`erase_if`, which implies in-place mutation matching C++20's `std::erase_if`. The
function has 6 active call sites in `HalfEdgeMesh.hpp` (lines 800, 813, 1094,
1096, 1117, 1119), all of which correctly capture the return value.

## Acceptance Criteria
- [ ] Either rename `detail::erase_if` to a name that reflects copy semantics
  (e.g. `detail::filter` or `detail::filtered_copy`), or replace all 6 call
  sites with inline `std::copy_if` / `std::ranges::filter_view` and remove the
  helper
- [ ] All tests pass
- [ ] No change in behavior at any call site

## Dependencies
None.
