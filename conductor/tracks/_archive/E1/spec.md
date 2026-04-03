# E1 — Flattening benchmark example

## GitHub Issue
https://github.com/educelab/OpenABF/issues/45

## Summary

Add a `BenchmarkFlattening` example app that measures wall-clock runtime for
multiple flattening configurations on one or more user-supplied meshes or
built-in wavy surfaces. Motivated by volume-cartographer#123 which benchmarked
ABF++ + LSCM (SparseLU vs CG) but noted HLSCM as a future alternative.

## Related
- Downstream: https://github.com/educelab/volume-cartographer/pull/123
- F1 (HLSCM): https://github.com/educelab/OpenABF/issues/5

## Acceptance Criteria
- [x] New `openabf_example_benchmark` target added to `examples/CMakeLists.txt`
- [x] Accepts mesh files as CLI arguments and/or `--builtin` for synthetic meshes
- [x] Times ABF++ angle optimization separately from each LSCM step
- [x] Benchmarks five configurations:
  - ABF++ + AngleBasedLSCM (SparseLU)
  - ABF++ + AngleBasedLSCM (CG `Lower|Upper`)
  - ABF++ + HierarchicalLSCM (LSCG)
  - ABF++ + HierarchicalLSCM (CG `Lower|Upper`)
- [x] Multi-thread columns via `--threads` option (explicit or auto-detect)
- [x] Prints results as a markdown table (faces, ABF time, per-solver columns)
- [x] Times in seconds with 2 decimal places
- [x] OOM and SolverException handled gracefully (prints "N/A")
- [x] Optional `--output-dir` writes flattened OBJ per algorithm

## Dependencies
- F1 (HLSCM): https://github.com/educelab/OpenABF/issues/5
