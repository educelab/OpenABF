# F2 Implementation Plan

## Phase 1 — Tests
1. Write tests with 2–3 pre-parameterized meshes
2. Assert packed charts do not overlap
3. Assert all UV coordinates are within [0,1]²

## Phase 2 — Design
1. Define the chart bounding box computation (min/max UV per mesh)
2. Choose packing algorithm: shelf-packing is simple and good enough for a first version
3. Decide on in-place UV update vs. separate UV buffer output

## Phase 3 — Implementation
1. Create `include/OpenABF/ChartPacking.hpp` (or add to a `Utils.hpp`)
2. Implement `PackCharts<MeshType>(std::vector<MeshType::Pointer>&)`
3. Add to `include/OpenABF/OpenABF.hpp`
4. Update `single_include.json`

## Phase 4 — Verify
- Run `ctest`
- Test with output from F3's `ExtractConnectedComponents` end-to-end
