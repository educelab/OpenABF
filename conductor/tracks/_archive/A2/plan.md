# A2 Implementation Plan

Implemented as part of P1. See P1/plan.md Phase 2 for the range type design.

## Specific steps
1. Add `auto edges()` and `auto edges() const` methods to `Face` that return the same
   iterable as `begin()`/`end()` but via a named range object
2. Update documentation examples to prefer `face->edges()`
3. Update internal usage in ABF, ABFPlusPlus, AngleBasedLSCM

## Verify
- Run `ctest`
