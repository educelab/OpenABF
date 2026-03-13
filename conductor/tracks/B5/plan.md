# B5 Implementation Plan

## Phase 1 — Test
Add a test case in `TestMeshIO.cpp` (or `TestMeshIOUtils.cpp`) that calls
`io_formats::is_file_type<io_formats::OBJ>` with a path that has no extension
and asserts the result is `false`.

## Phase 2 — Fix
Add an empty-string guard before `ext[0]`:
```cpp
if (ext.empty()) { return false; }
if (ext[0] == '.') { ext = ext.substr(1); }
```

## Phase 3 — Verify
- Run `ctest`
