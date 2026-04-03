# P4 Plan — IncompleteCholesky preconditioner for AngleBasedLSCM

## Phase 1: Benchmark

Add `CG(Diagonal)` and `CG(IC)` columns to `BenchmarkFlattening` and run
against real meshes to establish a three-way comparison vs. the current
`SparseLU` default.

- [ ] Add `AngleBasedLSCM` with `ConjugateGradient<Lower|Upper>` (Diagonal
  preconditioner) column to `BenchmarkFlattening.cpp`
- [ ] Add `AngleBasedLSCM` with `ConjugateGradient<Lower|Upper,
  IncompleteCholesky>` column to `BenchmarkFlattening.cpp`
- [ ] Run benchmark on real scroll meshes at multiple sizes
- [ ] Record results and determine if IC-CG beats SparseLU

## Phase 2: Default Change (conditional on Phase 1)

Only proceed if Phase 1 benchmarks show IC-CG is faster than SparseLU.

- [ ] Change `AngleBasedLSCM` default `Solver` template parameter from
  `SparseLU` to `ConjugateGradient<SparseMatrix<T>, Lower|Upper,
  IncompleteCholesky<T>>`
- [ ] Update `@tparam Solver` doc: explain IC vs Diagonal tradeoff, note
  SparseLU as reliable fallback for small meshes or ill-conditioned systems
- [ ] Run `ctest` — all tests must pass
- [ ] Run `git clang-format`
- [ ] Regenerate amalgamated header

## Phase 3: Conductor & GitHub

- [ ] Commit and push
- [ ] Update PR / close issue #47
- [ ] Update tracks.md
