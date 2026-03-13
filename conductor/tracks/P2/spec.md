# P2 — `std::map` used for interior-vertex index lookup in solver hot path

## GitHub Issue
https://github.com/educelab/OpenABF/issues/13

## Summary
`ABF::Compute` and `ABFPlusPlus::Compute` rebuild a `std::map<std::size_t, std::size_t>`
(`vIdx2vIntIdx`) every solver iteration to map vertex indices to interior-vertex
indices. Map lookup is O(log n); a flat vector with direct indexing would be O(1).

## Files
`include/OpenABF/ABF.hpp`, `include/OpenABF/ABFPlusPlus.hpp`

## Acceptance Criteria
- [ ] `vIdx2vIntIdx` replaced with a `std::vector<std::size_t>` of size `num_vertices()`
  indexed directly by vertex idx, or a dense array pre-built once per `Compute` call
- [ ] All existing tests pass

## Dependencies
- P1: if P1 redesigns the iteration data structures, the interior-vertex permutation
  may be naturally incorporated; P2 may be resolved as part of P1
