# P4 — IncompleteCholesky preconditioner for AngleBasedLSCM

## GitHub Issue
https://github.com/educelab/OpenABF/issues/47

## Summary

Evaluate `IncompleteCholesky`-preconditioned `ConjugateGradient` as a replacement
for the current `SparseLU` default solver in `AngleBasedLSCM`. If benchmarks show
IC-CG is faster than SparseLU (the current default) and faster than Diagonal-CG
(Eigen's CG default), change the `AngleBasedLSCM` default template parameter and
document the preconditioner choice as best practice in the `@tparam Solver` docs.

## Background

`AngleBasedLSCM` currently defaults to `SparseLU`, a direct solver. CG with
Eigen's default `DiagonalPreconditioner` (Jacobi) is already available by template
parameter but not the default. `IncompleteCholesky` is a better-fit preconditioner
for the LSCM normal equations (AᵀA is SPD), and is expected to reduce iteration
count by 3–5× vs. Diagonal on high-curvature meshes, potentially making CG+IC
faster than SparseLU for large meshes.

## Acceptance Criteria

- [ ] Benchmark compares three solver configurations on real meshes:
  - `SparseLU` (current default)
  - `ConjugateGradient` with `DiagonalPreconditioner` (Eigen CG default)
  - `ConjugateGradient` with `IncompleteCholesky`
- [ ] If IC-CG is faster than SparseLU: change `AngleBasedLSCM` default to IC-CG
- [ ] `@tparam Solver` doc updated to explain IC vs Diagonal tradeoff and recommend
  IC for large/high-curvature meshes
- [ ] All existing tests pass with updated default (if changed)

## Out of Scope

- Convenience type aliases
- Changes to HLSCM (warm-start already acts as effective preconditioner)
- Changes to LSCG path (IC only applies to square SPD systems)

## Dependencies

- E1 (BenchmarkFlattening): benchmark tool used for evaluation
