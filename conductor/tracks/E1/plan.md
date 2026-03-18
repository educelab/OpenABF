# E1 Plan — Flattening benchmark example

## Phase 1: Implement BenchmarkFlattening example
- [x] Create `examples/src/BenchmarkFlattening.cpp`
  - Parse mesh file paths and `--builtin`, `--threads`, `--output-dir` options
  - For each mesh (file or built-in wavy surface):
    - Run ABF++ once and record angle-optimization time
    - Time AngleBasedLSCM with SparseLU
    - Time AngleBasedLSCM with CG (`Lower|Upper`) at each thread count
    - Time HierarchicalLSCM with LSCG at each thread count
    - Time HierarchicalLSCM with CG (`Lower|Upper`) at each thread count
  - Print markdown table with all columns
  - OOM/SolverException handling (prints "N/A", continues)
  - Optional `--output-dir` writes flattened OBJ per algorithm
- [x] Add `openabf_example_benchmark` to `examples/CMakeLists.txt`
- [x] Regenerate amalgamated header (no header changes needed — example only)

## Phase 2: Conductor & GitHub
- [x] Commit and push (landed on `feat/f1-hlscm` branch, PR #44)
- [x] Update tracks.md
