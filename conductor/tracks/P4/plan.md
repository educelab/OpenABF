# P4 Implementation Plan

## Phase 1: Benchmark
Add `CG(Diagonal)` and `CG(IC)` columns to `BenchmarkFlattening` and run against real meshes to establish a three-way comparison vs. the current `SparseLU` default.

- [ ] 1.1 Add `AngleBasedLSCM` with `ConjugateGradient<Lower|Upper>` (Diagonal preconditioner) column to `BenchmarkFlattening.cpp`
- [ ] 1.2 Add `AngleBasedLSCM` with `ConjugateGradient<Lower|Upper, IncompleteCholesky>` column to `BenchmarkFlattening.cpp`
- [ ] 1.3 Run benchmark on real scroll meshes at multiple sizes
- [ ] 1.4 Record results and determine if IC-CG beats SparseLU

## Phase 2: Default Change (conditional on Phase 1)
Only proceed if Phase 1 benchmarks show IC-CG is faster than SparseLU.

- [ ] 2.1 Change `AngleBasedLSCM` default `Solver` template parameter from `SparseLU` to `ConjugateGradient<SparseMatrix<T>, Lower|Upper, IncompleteCholesky<T>>`
- [ ] 2.2 Update `@tparam Solver` doc: explain IC vs Diagonal tradeoff, note SparseLU as reliable fallback for small meshes or ill-conditioned systems
- [ ] 2.3 Run `ctest` — all tests must pass
- [ ] 2.4 Run `git clang-format`
- [ ] 2.5 Regenerate amalgamated header

## Phase 3: Conductor & GitHub
- [ ] 3.1 Commit and push
- [ ] 3.2 Update PR / close issue #47
- [ ] 3.3 Update tracks.md
