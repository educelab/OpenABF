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
- [x] 2.1 Synthetic 2D charts: bbox computation correctness (min/max per mesh)
- [x] 2.2 Assert packed chart bounding boxes do not overlap (padding respected)
- [x] 2.3 Assert absolute-mode preserves relative chart sizes (no per-chart distortion)
- [x] 2.4 Assert normalize=true fits all UVs within [0,1]² via single global scale
- [x] 2.5 Assert returned PackResult extent bounds all packed charts
- [x] 2.6 Degenerate cases: empty list, single chart, zero-area chart, null/empty throw
- [x] 2.7 End-to-end: tear → extract_connected_components → LSCM → PackCharts, and
          verify per-wedge recovery via (face_map[f], vertex_map[corner.vertex.idx])

## Phase 3: Implementation
- [x] 3.1 Create `include/OpenABF/ChartPacking.hpp` with PackOptions, PackResult
- [x] 3.2 Implement per-chart bbox + sqrt-area target width + shelf placement
- [x] 3.3 Implement absolute (translate-only) and normalize (global uniform scale) modes
- [x] 3.4 Implement padding, degenerate-input handling, static_assert(Dim>=2)
- [x] 3.5 Document the vertex-identity per-wedge recipe in the header + complexity notes
- [x] 3.6 Add include to `include/OpenABF/OpenABF.hpp`
- [x] 3.7 Update single-header via amalgamation script (single_include.json unchanged —
          it already tracks OpenABF.hpp transitively)

## Phase 4: Verify
- [x] 4.1 Run `ctest` — all suites pass (incl. OpenABF_TestChartPacking, OpenABF_TestMeshMerge)
- [x] 4.2 Run clang-format on changed files
- [x] 4.3 Confirm single-header build compiles and runs

## Phase 5: MergeMeshes helper (added during review)
Rationale: the inline atlas merge in the example severs the back-map chain.
MergeMeshes is the inverse of extract_connected_components — it returns
provenance maps so merged → chart → M' composition keeps working.
- [x] 5.1 Tests first: concatenation counts, vertex/face provenance, null/empty
          throw, round-trip extract→merge recovers original (torn-mesh) identity
- [x] 5.2 Implement `MergeMeshes<MeshType>` → `MergedMesh{mesh, vertex_source, face_source}`
          in `include/OpenABF/MeshMerge.hpp` (preserves vertex traits/positions;
          edge/face traits default-constructed)
- [x] 5.3 Wire into OpenABF.hpp + multiheader install list; regenerate single header
- [x] 5.4 Switch MultiChartFlatten example to use MergeMeshes
- [x] 5.5 Verify: full ctest, single-header build, install-test all pass

## Phase 6: Perimeter padding (added during review)
Rationale: review question "shouldn't we add padding around the packed
charts?". The original layout applied `padding` only as a gutter *between*
charts — perimeter charts still touched the atlas boundary (left/bottom at the
origin, rightmost/topmost at the extent). For a texture atlas this lets edge
charts bleed across the boundary/seam under filtering, mipmapping, or wrap
addressing. Resolution (user-confirmed): inset the whole layout so `padding`
surrounds every chart on all four sides; keep the library default `padding = 0`
and instead set a visible padding in the example.
- [x] 6.1 Tests first: assert padding insets charts from the atlas perimeter
          (new PaddingSurroundsChartsAtPerimeter); update single-row extent
          expectation (pad + w0 + pad + w1 + pad)
- [x] 6.2 Implement perimeter inset: cursor starts/wraps at `pad`; add `pad` to
          far extents; normalize fits the padded atlas into [0,1]²
- [x] 6.3 Update header docs (padding surrounds charts; atlas lower corner stays
          at origin) + spec Decision 5 / acceptance criteria
- [x] 6.4 Set a visible `padding` in the MultiChartFlatten example
- [x] 6.5 Regenerate single header; verify full ctest, example run, single-header
          build, clang-format all pass
