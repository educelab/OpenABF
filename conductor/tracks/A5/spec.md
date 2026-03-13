# A5 — No configurable pinned edge selection for LSCM

## GitHub Issue
https://github.com/educelab/OpenABF/issues/2

## Summary
`AngleBasedLSCM` always selects the first boundary edge as the pinned edge, and maps
it to the closest XY axis. Better heuristics (e.g., the longest boundary edge, a
user-specified edge, or Blender's balanced approach) can significantly improve
parameterization quality by reducing distortion.

## File
`include/OpenABF/AngleBasedLSCM.hpp`

## Related
GitHub issue #2

## Acceptance Criteria
- [ ] The pinned vertex/edge can be overridden by the caller
- [ ] The default behavior (first boundary edge) is preserved when no override is given
- [ ] An optional heuristic is provided (e.g. longest boundary edge)
- [ ] Tests cover at least the default and one override strategy
- [ ] All existing tests pass

## Proposed API
A strategy object or optional vertex-index pair:
```cpp
// Default: auto-select
LSCM::Compute(mesh);

// Caller specifies pinned vertices by index
LSCM::Compute(mesh, /*p0=*/0, /*p1=*/3);
```

## Dependencies
None.
