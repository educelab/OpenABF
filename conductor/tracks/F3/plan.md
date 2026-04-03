# F3 Implementation Plan

## Phase 1 — Tests
1. Test `ExtractConnectedComponents` on:
   - A single-CC mesh (returns one element, same geometry)
   - A two-CC mesh (returns two elements, correct vertex maps)
   - A mesh with pre-computed angle traits (traits are preserved in extracted meshes)
2. Test `ParameterizeConnectedComponents` runs without error on a two-CC mesh

## Phase 2 — ExtractConnectedComponents
1. Create `include/OpenABF/HalfEdgeMeshUtils.hpp`
2. Implement `ExtractConnectedComponents`:
   - Call `mesh->connected_components()` to get face groups
   - For each group, collect unique vertex indices, build a local index map
   - Create a new `MeshType` mesh, insert vertices (copying traits via copy ctor)
   - Insert faces using re-mapped indices (copying edge traits)
   - Call `update_boundary()` on extracted mesh
   - Build `original_idx` back-map vector
3. Return the vector of pairs

## Phase 3 — ParameterizeConnectedComponents
1. Implement the convenience wrapper:
   - Call `ExtractConnectedComponents`
   - For each `(extracted_mesh, back_map)`:
     - Run `AngleOptimizer::Compute(extracted_mesh)` if optimizer is not `void`
     - Run `Parameterizer::Compute(extracted_mesh)`
     - Write UV positions back to original mesh using `back_map`
2. Add to `include/OpenABF/OpenABF.hpp` and `single_include.json`

## Phase 4 — Verify
- Run `ctest`
- End-to-end test: split_path → extract → parameterize → pack (with F2)
