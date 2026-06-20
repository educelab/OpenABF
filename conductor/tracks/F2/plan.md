# F2 Implementation Plan

Design resolved 2026-06-20 (see spec.md → Design Decisions).

## Phase 1: Design freeze
- [x] 1.1 Resolve scope: geometry-only, in-place pos mutation (no UVMap class)
- [x] 1.2 Resolve per-wedge recipe: key on vertex identity via back-maps, not corner position
- [x] 1.3 Resolve scaling: absolute by default, opt-in global-uniform normalize to [0,1]²
- [x] 1.4 Resolve layout: shelf packing, ~square sqrt-area target width (overridable)
- [x] 1.5 Resolve API: PackOptions{normalize,target_width,padding}, PackResult{min,max}
- [x] 1.6 Resolve edge cases: empty list no-op, zero-area placed, null/empty throws, Dim>=2

## Phase 2: Tests (write first)
- [ ] 2.1 Synthetic 2D charts: bbox computation correctness (min/max per mesh)
- [ ] 2.2 Assert packed chart bounding boxes do not overlap (padding respected)
- [ ] 2.3 Assert absolute-mode preserves relative chart sizes (no per-chart distortion)
- [ ] 2.4 Assert normalize=true fits all UVs within [0,1]² via single global scale
- [ ] 2.5 Assert returned PackResult extent bounds all packed charts
- [ ] 2.6 Degenerate cases: empty list, single chart, zero-area chart, null/empty throw
- [ ] 2.7 End-to-end: tear → extract_connected_components → LSCM → PackCharts, and
          verify per-wedge recovery via (face_map[f], vertex_map[corner.vertex.idx])

## Phase 3: Implementation
- [ ] 3.1 Create `include/OpenABF/ChartPacking.hpp` with PackOptions, PackResult
- [ ] 3.2 Implement per-chart bbox + sqrt-area target width + shelf placement
- [ ] 3.3 Implement absolute (translate-only) and normalize (global uniform scale) modes
- [ ] 3.4 Implement padding, degenerate-input handling, static_assert(Dim>=2)
- [ ] 3.5 Document the vertex-identity per-wedge recipe in the header + complexity notes
- [ ] 3.6 Add include to `include/OpenABF/OpenABF.hpp`
- [ ] 3.7 Update `single_include.json` and run amalgamation script

## Phase 4: Verify
- [ ] 4.1 Run `ctest` — all tests pass
- [ ] 4.2 Run clang-format on changed files
- [ ] 4.3 Confirm single-header build matches multi-header
