# P7 Implementation Plan

Strict TDD (per `conductor/workflow.md`): failing tests first, then the
implementation, then docs/amalgamation, with manual verification at each phase
checkpoint.

**Blocked until PR #99 (F2) merges.** Phase 0 is the gate.

## Phase 0: Gate on F2
- [ ] 0.1 Confirm PR #99 is merged and `include/OpenABF/ChartPacking.hpp` plus
      `tests/src/TestChartPacking.cpp` are on `develop`.
- [ ] 0.2 Branch `p7-packing-square-wrap` from updated `develop`; open a draft PR
      linked to issue #104 via GitHub's Development field.

### Verification
- [ ] `ctest` is green on the fresh branch before any P7 change (clean baseline).

## Phase 1: Pin current behavior with failing tests
Write the tests that the current heuristic fails, and record what it actually
produces so the improvement is measurable.

- [ ] 1.1 Add `DefaultWrapPacksTallChartPairSideBySide` to
      `tests/src/TestChartPacking.cpp`: two 1 × 2 charts,
      `minimize_bounding_box = false`, `padding = 0.1`, no `target_width`; expect
      extent 2.3 × 2.2. (Currently 1.2 × 4.3.)
- [ ] 1.2 Add `DefaultWrapKeepsAtlasRoughlySquare`: chart sets where near-square is
      achievable (4, 6, and 9 unit charts; `padding = 0` and `padding = 0.5`);
      assert `max(W, H) / min(W, H)` is within a tight bound. Excludes counts
      (e.g. 2 identical squares) where no layout beats aspect 2.
- [ ] 1.3 Add `DefaultWrapPaddingDoesNotForceEarlyWrap`: for a set of identical
      charts, the charts-per-shelf count with `padding > 0` matches the
      `padding = 0` count.
- [ ] 1.4 Add `ExplicitTargetWidthWrapsAtGivenWidth`: pins the override path's wrap
      case (the existing `PaddingSeparatesChartsInSingleRow` only pins the
      no-wrap case), and `OversizedChartPlacedAloneOnShelf`: a chart wider than
      `target_width` is placed at a shelf start rather than wrapping an empty
      shelf.
- [ ] 1.5 Run `ctest`; record which new tests fail and the extents the current
      heuristic produces, in this plan.

### Verification
- [ ] 1.1–1.3 fail against unmodified `PackCharts`; 1.4 passes (override path is
      unchanged by P7 and must stay that way).
- [ ] Measured pre-change extents recorded below.

### Phase 1 results
_(to be filled in: per-case current vs. expected extents)_

## Phase 2: Squareness-based wrap
- [ ] 2.1 In `PackCharts` (`include/OpenABF/ChartPacking.hpp`), add a local
      lexicographic cost on `(max(W, H), W * H)` and compute the stay/wrap
      candidate extents from `cursorX`, `cursorY`, `shelfHeight`, `atlasMaxX`,
      `atlasMaxY` as in the spec.
- [ ] 2.2 Replace the default wrap test with `cursorX > pad and cost(stay) >
      cost(wrap)`; keep the `cursorX > pad` empty-shelf guard so oversized charts
      are still placed alone. Strict `>` so ties keep the chart on the shelf.
- [ ] 2.3 Keep `opts.target_width` as an outright override on its existing
      `cursorX + width[i] > targetWidth` test; delete the now-dead
      `sqrt(Σ (w+pad)(h+pad))` `areaSum` loop from the default path.
- [ ] 2.4 Run `ctest` — Phase 1 tests 1.1–1.3 now pass, 1.4 still passes, and all
      pre-existing `ChartPacking` tests pass **unmodified**.
- [ ] 2.5 Confirm `PaddingSurroundsChartsAtPerimeter` still spans ≥ 2 shelves under
      the new decision (its stated purpose); if its 5 charts now fit one row,
      adjust the chart set — not the assertions — so the top/right perimeter
      margins stay covered, and update its comment.

### Verification
- [ ] Full `ctest` green; no existing test assertion weakened or deleted.
- [ ] Non-overlap and perimeter-inset invariants still hold (`ChartsDoNotOverlap`,
      `ExtentBoundsAllCharts`, `PaddingSurroundsChartsAtPerimeter`).

## Phase 3: Documentation and single header
- [ ] 3.1 Rewrite the `PackOptions::target_width` Doxygen: no `sqrt(total area)`
      default; describe the squareness-based wrap when unset, and the override's
      exact semantics (unpadded right edge vs. the given width).
- [ ] 3.2 Review the `PackCharts` shelf-layout comment block and `@par Complexity`
      — describe the new decision and confirm `O(n log n)` still holds.
- [ ] 3.3 `git clang-format` on changed files, re-stage.
- [ ] 3.4 Regenerate the amalgamated header:
      `python3 thirdparty/amalgamate/amalgamate.py -c single_include.json -s .`
      (no new headers introduced, so no `single_include.json` edit expected).
- [ ] 3.5 Build docs and confirm no new Doxygen warnings.

### Verification
- [ ] `ctest` green against the regenerated single header (tests include only
      `OpenABF.hpp`).
- [ ] `git diff` on `single_include/` reflects exactly the `ChartPacking.hpp`
      change.

## Phase 4: Example verification and integration
- [ ] 4.1 Build and run `examples/src/MultiChartFlatten.cpp`; confirm the reported
      normalized extent is near-square instead of `[0, 0] -> [0.27907, 1]`, and
      record the before/after values in this plan.
- [ ] 4.2 Confirm all other examples still build and run.
- [ ] 4.3 Commit with Conventional Commits (`test:` for Phase 1, `fix:` or `perf:`
      for Phase 2, `docs:`/`chore:` for Phase 3); push.
- [ ] 4.4 Mark the PR ready, update its title and description with the full diff
      summary, and confirm it closes #104.
- [ ] 4.5 Update `conductor/tracks.md`: P7 → closed/archived with the PR number,
      and refresh the `Updated` date.

### Phase 4 results
_(to be filled in: MultiChartFlatten extent before/after)_

## Final Verification
- [ ] All spec acceptance criteria met.
- [ ] `ctest` fully green; new tests cover the default heuristic (the gap #104
      identified), the padding-independence of shelf occupancy, and both
      `target_width` override paths.
- [ ] `git clang-format` clean; single header regenerated and committed.
- [ ] Docs updated with no new warnings; example output verified.
- [ ] PR reviewed and ready to merge.
