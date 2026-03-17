# E1 — Flattening benchmark example

## GitHub Issue
https://github.com/educelab/OpenABF/issues/45

## Summary

Add a `BenchmarkFlattening` example app that measures wall-clock runtime for
multiple flattening configurations on one or more user-supplied meshes. Motivated
by volume-cartographer#123 which benchmarked ABF++ + LSCM (SparseLU vs CG) but
noted HLSCM as a future alternative.

## Related
- Downstream: https://github.com/educelab/volume-cartographer/pull/123
- F1 (HLSCM): https://github.com/educelab/OpenABF/issues/5

## Acceptance Criteria
- [ ] New `openabf_example_benchmark` target added to `examples/CMakeLists.txt`
- [ ] Accepts one or more mesh files as CLI arguments
- [ ] Times ABF++ angle optimization separately from each LSCM step
- [ ] Benchmarks three configurations:
  - ABF++ + AngleBasedLSCM (SparseLU)
  - ABF++ + AngleBasedLSCM (LSCG)
  - ABF++ + HierarchicalLSCM (LSCG)
- [ ] Prints results as a markdown table (faces, ABF time, LSCM-SparseLU, LSCM-LSCG, HLSCM)
- [ ] Times in seconds with 2 decimal places
