# P6 — Precompute sin/cos values before HLSCM face assembly loop

## GitHub Issue
https://github.com/educelab/OpenABF/issues/66

## Status
**Blocked on [A8](../A8/spec.md) (#64).**

## Re-evaluate after A8
Before starting work on this track, confirm it is still worth doing:
- A8 extracts the HLSCM face-assembly loop into `detail::lscm::buildSystem`.
  The shape of that loop after A8 may already eliminate the redundancy P6
  targets (e.g. if angles are walked once per face rather than per edge).
- Re-run the benchmark on a representative mesh against the post-A8 baseline.
  If the sin/cos cost is no longer a measurable hot spot, close P6 as obsolete
  rather than implementing it.
- If P6 is still warranted, scope it against the post-A8 loop, not the
  pre-A8 version — the file, function, and surrounding code will have moved.

## Dependencies
- **Hard:** [A8](../A8/spec.md) (#64) must land first. A8 extracts the LSCM
  system-building logic (including the face-assembly loop P6 modifies) into
  `detail::lscm::buildSystem`. Doing P6 first would force a rewrite during A8's
  merge.

## Summary
Precompute trigonometric sin/cos values used in the HLSCM face assembly loop to avoid redundant computation per iteration.

## Acceptance Criteria
- [ ] Identify all sin/cos calls inside HLSCM face assembly loop
- [ ] Precompute values before the loop and store in local arrays
- [ ] Verify no regression in parameterization correctness
- [ ] Benchmark to confirm performance improvement
