# E1 Plan — Flattening benchmark example

## Phase 1: Implement BenchmarkFlattening example
- [ ] Create `examples/src/BenchmarkFlattening.cpp`
  - Parse mesh file paths from argv
  - For each mesh file:
    - Load mesh, print face/vertex count
    - Run ABF++ once and record angle-optimization time
    - Time AngleBasedLSCM with SparseLU (re-clone mesh between runs)
    - Time AngleBasedLSCM with LSCG
    - Time HierarchicalLSCM with LSCG
  - Print markdown table row per mesh
- [ ] Add `openabf_example_benchmark` to `examples/CMakeLists.txt`
- [ ] Regenerate amalgamated header (no header changes needed — example only)

## Phase 2: Conductor & GitHub
- [ ] Commit and push on `feat/e1-benchmark` branch
- [ ] Open PR linking issue #45
- [ ] Update tracks.md
