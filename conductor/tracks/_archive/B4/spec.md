# B4 — `Vec` binary `*` and `/` operators inconsistent with `*=` and `/=`

## GitHub Issue
https://github.com/educelab/OpenABF/issues/9

## Summary
`Vec::operator*=(T2)` and `operator/=(T2)` require an arithmetic scalar, but the
corresponding binary `operator*(Vec, Vector)` and `operator/(Vec, Vector)` accept
any `Vector` type and call `*=` / `/=` with it — which will fail to compile or
silently misbehave when `Vector` is not a scalar.

## File
`include/OpenABF/Vec.hpp`

## Acceptance Criteria
- [ ] Binary `operator*` and `operator/` are constrained to arithmetic scalar types,
  matching `operator*=` and `operator/=`
- [ ] Tests for scalar multiplication/division pass
- [ ] No unintended element-wise vector multiplication is introduced (this is not a
  planned operation for `Vec`)

## Dependencies
None.
