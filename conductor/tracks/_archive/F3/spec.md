# F3 — Multi-component extraction and parameterization pipeline

## GitHub Issue
https://github.com/educelab/OpenABF/issues/19

## Summary
There is no way to feed a mesh with multiple connected components (e.g. after
`split_path` calls) directly into the parameterization pipeline. ABF and LSCM both
assume a single connected component. A two-layer API is proposed:

1. **Low-level extractor** — `ExtractConnectedComponents` returns a vector of
   `(extracted_mesh_ptr, original_vertex_index_map)` pairs
2. **High-level convenience helper** — `ParameterizeConnectedComponents` runs
   the full angle-optimizer + parameterizer pipeline on each CC and writes UV
   results back to the original mesh

## Acceptance Criteria
- [ ] `ExtractConnectedComponents<MeshType>(mesh)` returns
  `std::vector<std::pair<MeshType::Pointer, std::vector<std::size_t>>>`
  where `pair.second[extracted_idx] == original_idx`
- [ ] Extracted meshes are fully independent deep copies (including edge/vertex traits)
- [ ] Vertex indices in extracted meshes are re-densified (0..N-1)
- [ ] `ParameterizeConnectedComponents<AngleOptimizer, Parameterizer>(mesh)` runs the
  full pipeline and writes UV coordinates back to the original mesh
- [ ] A single-component mesh passes through unchanged
- [ ] Tests cover: single CC (passthrough), two CCs, trait preservation
- [ ] Lives in a new header (e.g. `HalfEdgeMeshUtils.hpp`) not in `HalfEdgeMesh.hpp`
- [ ] Included in `OpenABF.hpp`

## Dependencies
- P1 (completing the iteration redesign first is recommended for performance;
  extraction uses iteration heavily)
- F2 (F3's output feeds into F2's chart packing)
