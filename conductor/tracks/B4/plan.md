# B4 Implementation Plan

## Phase 1 — Tests
1. In `tests/src/TestVec.cpp`, add tests confirming scalar `*` and `/` work correctly
2. Confirm the test fails to compile (or produces wrong results) before the fix

## Phase 2 — Fix
Change the binary `operator*` and `operator/` friend signatures to require an
arithmetic scalar, matching `operator*=`:

```cpp
template <typename T2, std::enable_if_t<std::is_arithmetic<T2>::value, bool> = true>
friend Vec operator*(Vec lhs, const T2& rhs) { lhs *= rhs; return lhs; }

template <typename T2, std::enable_if_t<std::is_arithmetic<T2>::value, bool> = true>
friend Vec operator/(Vec lhs, const T2& rhs) { lhs /= rhs; return lhs; }
```

## Phase 3 — Verify
- Run `ctest`
- Confirm that passing a `Vec` as the RHS of `*` or `/` now fails to compile
  (add a `static_assert` test if feasible)
