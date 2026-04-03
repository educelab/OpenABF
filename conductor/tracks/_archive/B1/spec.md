# B1 — ABF lambda update index wrong for meshes with 2+ interior vertices

## GitHub Issue
https://github.com/educelab/OpenABF/issues/6

## Summary
The interior-vertex lambda update loop in `ABF::Compute` uses `idx + intIdx` where
both variables increment each iteration, producing wrong offsets into the delta vector.
Only affects meshes with 2+ interior vertices; the single-interior-vertex pyramid test
masks this completely.

## File
`include/OpenABF/ABF.hpp` (~line 398)

## Reproduction
Any mesh with ≥ 2 interior vertices passed through `ABF::Compute` will produce
incorrect Lagrange multiplier updates and thus a wrong parameterization.

## Acceptance Criteria
- [ ] Tests exist for ABF on a mesh with ≥ 3 interior vertices (see T1)
- [ ] `v->lambda_plan` and `v->lambda_len` are updated using the correct delta indices
- [ ] All existing parameterization tests continue to pass
- [ ] ABF produces correct results on a flat grid mesh (verifiable analytically)

## Dependencies
- T1 must be completed first (need the multi-interior-vertex test mesh)

## Notes
The fix is to replace the running `idx + intIdx` pattern with a fixed base offset:
```cpp
auto base = edgeCnt + faceCnt;
// then use base + intIdx and base + vIntCnt + intIdx
```
`ABFPlusPlus` gets this right using `faceCnt + intIdx` directly — use it as reference.
