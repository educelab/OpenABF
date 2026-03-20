# T3 Implementation Plan

## Phase 1 — Tests
In `tests/src/TestMeshIO.cpp`:

1. `TEST(MeshIO, OBJ_RoundTrip)`:
   - Construct the pyramid mesh
   - Write to `std::ostringstream` using `io_formats::OBJ::Write`
   - Read back from `std::istringstream` using `io_formats::OBJ::Read`
   - Assert vertex count, face count, and all vertex positions match

2. `TEST(MeshIO, PLY_RoundTrip)`:
   - Same, using `io_formats::PLY::Write` and `PLY::Read`

## Phase 2 — Verify
- Run `ctest`
