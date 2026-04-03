# B5 — `is_file_type` accesses `ext[0]` without bounds check

## GitHub Issue
https://github.com/educelab/OpenABF/issues/10

## Summary
In `MeshIOFormats.hpp`, `is_file_type` calls `ext[0]` on the result of
`path.extension().string()`. If the path has no extension, `ext` is empty and
`ext[0]` is undefined behavior.

## File
`include/OpenABF/MeshIOFormats.hpp` (line 22)

## Acceptance Criteria
- [ ] A path with no extension (e.g. `"mesh"`) does not invoke UB
- [ ] `is_file_type` returns `false` for a path with no extension
- [ ] All existing IO tests pass

## Dependencies
None.
