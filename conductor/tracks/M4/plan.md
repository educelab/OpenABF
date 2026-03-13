# M4 Implementation Plan

## Phase 1 — Fix
Move the `operator<<` definition inside the `OpenABF` namespace, after the `Vec` class
definition and before the closing brace.

## Phase 2 — Verify
- Run `ctest`
- Confirm existing streaming uses in tests still compile
