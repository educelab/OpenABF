# B6 — PLY reader: `vmap` may silently read wrong column if x/y/z absent

## GitHub Issue
https://github.com/educelab/OpenABF/issues/11

## Summary
The PLY reader initializes `std::array<std::size_t, 3> vmap{}` to `{0, 0, 0}`.
If the PLY header is missing one or more of the `x`, `y`, or `z` vertex properties,
those indices stay at 0 and the reader silently reads the wrong data column.

## File
`include/OpenABF/MeshIOFormats.hpp`

## Acceptance Criteria
- [ ] After parsing the header, the reader validates that all three of `x`, `y`, `z`
  were found in the vertex element properties
- [ ] A PLY file missing any of these properties throws a `std::runtime_error`
- [ ] All existing PLY IO tests pass

## Dependencies
None.
