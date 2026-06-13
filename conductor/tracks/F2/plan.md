# F2 Implementation Plan

## Phase 1: Tests
- [ ] 1.1 Write tests with 2–3 pre-parameterized meshes
- [ ] 1.2 Assert packed charts do not overlap
- [ ] 1.3 Assert all UV coordinates are within [0,1]²

## Phase 2: Design
- [ ] 2.1 Define the chart bounding box computation (min/max UV per mesh)
- [ ] 2.2 Choose packing algorithm: shelf-packing is simple and good enough for a first version
- [ ] 2.3 Decide on in-place UV update vs. separate UV buffer output

## Phase 3: Implementation
- [ ] 3.1 Create `include/OpenABF/ChartPacking.hpp` (or add to a `Utils.hpp`)
- [ ] 3.2 Implement `PackCharts<MeshType>(std::vector<MeshType::Pointer>&)`
- [ ] 3.3 Add to `include/OpenABF/OpenABF.hpp`
- [ ] 3.4 Update `single_include.json`

## Phase 4: Verify
- [ ] 4.1 Run `ctest`
- [ ] 4.2 Test with output from F3's `ExtractConnectedComponents` end-to-end
