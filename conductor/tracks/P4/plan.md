# P4 Implementation Plan

## Phase 1: Benchmark
Add `CG(Diagonal)` and `CG(IC)` columns to `BenchmarkFlattening` and run against real meshes to establish a three-way comparison vs. the current `SparseLU` default.

- [x] 1.1 Add `AngleBasedLSCM` with `ConjugateGradient<Lower|Upper>` (Diagonal preconditioner) column to `BenchmarkFlattening.cpp` — the existing `LSCM CG` column was already this (Eigen's CG default preconditioner is `DiagonalPreconditioner`); renamed to `LSCM CG-Diag` for clarity.
- [x] 1.2 Add `AngleBasedLSCM` with `ConjugateGradient<Lower|Upper, IncompleteCholesky>` column to `BenchmarkFlattening.cpp` (`LSCM CG-IC`).
- [x] 1.3 Run benchmark on built-in wavy + sphere mesh sequences at 50k–200k faces.
- [x] 1.4 Record results — IC-CG does NOT beat SparseLU; in fact much slower.

### Phase 1 results (1 thread, no OpenMP, FloatT=double)

| Mesh        | Faces  | SparseLU | CG-Diag | CG-IC  | CG-IC / SparseLU |
|-------------|--------|----------|---------|--------|------------------|
| wavy~50k    | 49928  | 0.32s    | 3.80s   | 14.36s | 45×              |
| wavy~100k   | 99458  | 0.91s    | 13.89s  | 69.17s | 76×              |
| sphere~50k  | 50880  | 0.52s    | 2.24s   | 7.57s  | 15×              |

The gap *widens* with mesh size: IC-CG is being dominated by IncompleteCholesky's factorization setup cost and a high per-iteration apply cost on these well-conditioned LSCM normal equations.  IC-CG is also markedly slower than Diagonal-CG, so on iterative-solver workloads `Diagonal` remains the better default preconditioner.

## Phase 2: Default Change (conditional on Phase 1)
**Not triggered.** Phase 1 acceptance criterion ("If IC-CG is faster than SparseLU") was not met. `AngleBasedLSCM`'s default `Solver` stays `SparseLU`.

- [x] 2.2 Update `@tparam Solver` doc — recorded that SparseLU is the benchmarked-fastest default in the 50k–200k face range, and that IC is available but currently underperforms Diagonal on the LSCM normal equations.
- [x] 2.3 Run `ctest` — all 43 parameterization tests pass on `p4-ic-cg-preconditioner` against the modified BenchmarkFlattening + updated doc.
- [x] 2.4 Run `git clang-format` — no changes needed.
- [x] 2.5 Regenerate amalgamated header — single_include updated.

## Phase 3: Conductor & GitHub
- [x] 3.1 Commit and push.
- [x] 3.2 PR #91 opened and merged against #47; benchmark/conclusion/doc update recorded; issue #47 closed as "investigated, default unchanged".
- [x] 3.3 Update tracks.md (P4 moved to archived).
