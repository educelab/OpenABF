# P7 — PackCharts wastes atlas area: shelf-wrap test and sqrt-area target width disagree on padding

## GitHub Issue
https://github.com/educelab/OpenABF/issues/104

## Summary
`PackCharts`'s default shelf-width heuristic and its shelf-wrap test account for
`padding` differently, so charts wrap to a new shelf earlier than the heuristic
intends. The layout stays valid and non-overlapping — this is not a correctness
bug — but it wastes atlas area, and at low chart counts a single premature wrap
turns a would-be square atlas into a narrow column. Introduced with F2 (#18,
PR #99).

Fix by removing the estimate/measurement mismatch entirely: wrap when wrapping
leaves the running atlas extent closer to square than staying on the current
shelf, rather than when a row would exceed an estimated width. An explicit
`opts.target_width` continues to override with today's semantics.

## Problem description
The target width is computed from **padded** chart areas
(`include/OpenABF/ChartPacking.hpp:336-342`):

```cpp
areaSum += (width[i] + pad) * (height[i] + pad);
targetWidth = std::sqrt(areaSum);
```

but the wrap test measures the **unpadded** right edge of the chart being placed
against a cursor that already starts at the perimeter inset `pad`
(`include/OpenABF/ChartPacking.hpp:358`):

```cpp
if (cursorX > pad and cursorX + width[i] > targetWidth) { /* wrap */ }
```

`sqrt(Σ (w+pad)(h+pad))` under-counts padding: a row of *k* charts spends *k+1*
horizontal gutters (the perimeter inset at both ends plus the interior ones), not
*k*. The estimate lands just below the width the charts actually need and the row
wraps early.

Worked for two 1 × 2 charts with `padding = 0.1`:

| step | value |
| ---- | ----- |
| `areaSum` | `2 × (1.1 × 2.1)` = 4.62 |
| `targetWidth` | `sqrt(4.62)` ≈ 2.149 |
| chart A at `cursorX = 0.1` | `0.1 + 1 = 1.1 ≤ 2.149` ✔ |
| cursor advances | `0.1 + 1 + 0.1 = 1.2` |
| chart B | `1.2 + 1 = 2.2 > 2.149` → **wrap** (by 0.05) |

Side by side the pair needs `pad + 1 + pad + 1 + pad` = 2.3 wide by
`pad + 2 + pad` = 2.2 tall (area 5.06, more than the estimated 4.62), giving an
almost exactly square atlas. Instead `examples/src/MultiChartFlatten.cpp` (two
charts ≈ 1 × 2 after bounding-box minimization, `padding = 0.1`,
`normalize = true`) reports a single stacked column:

```
Packed atlas extent: [0, 0] -> [0.27907, 1]     // 1.2 x 4.3 before normalization
```

## Why this matters
The wasted area is charged directly to the consumer as soon as the packed extent
is fit to a square texture — which is exactly what `opts.normalize` does, scaling
by `1 / max(atlasW, atlasH)`. For the two-chart example above the atlas occupies
`1.2 × 4.3 = 5.16` units inside a `4.3 × 4.3 = 18.5` square, a 28% fill; the
side-by-side layout fills `5.06 / 5.29` = 96% of its square. The effect is worst
at 2–4 charts and shrinks as chart count grows and rows fill.

## Approach — wrap on squareness
Replace the fixed-width wrap test in the default path with a greedy comparison of
the two candidate atlas extents at each placement, and drop the now-unused
sqrt-area estimate. At the decision point for chart `i` the loop already holds
`cursorX`, `cursorY`, `shelfHeight`, `atlasMaxX`, and `atlasMaxY`, so both
candidates are O(1):

```cpp
// stay on this shelf
const T stayW = std::max(atlasMaxX, cursorX + width[i]) + pad;
const T stayH = std::max(atlasMaxY, cursorY + height[i]) + pad;
// wrap to a new shelf
const T wrapY = cursorY + shelfHeight + pad;
const T wrapW = std::max(atlasMaxX, pad + width[i]) + pad;
const T wrapH = std::max(atlasMaxY, wrapY + height[i]) + pad;
if (cursorX > pad and cost(stayW, stayH) > cost(wrapW, wrapH)) { /* wrap */ }
```

`cost` is lexicographic on `(max(W, H), W * H)`: minimize the side of the
enclosing square first — the quantity `normalize` divides by and the one a square
texture pays for — then break ties toward the smaller total extent. Strict `>`
means an exact tie keeps the chart on the current shelf, so the layout prefers
fewer shelves and stays deterministic for identical inputs.

The existing `cursorX > pad` guard is retained, so an empty shelf never wraps and
a chart wider than any sensible row is still placed alone at a shelf start.
Charts are still visited in descending-height order, so `height[i] <= shelfHeight`
for every non-first chart on a shelf and `stayH` reduces to the current
`atlasMaxY + pad` in practice; the `std::max` is kept for clarity and safety.

Expected extents (`minimize_bounding_box = false`, `normalize = false`):

| charts | padding | current | after P7 |
| ------ | ------- | ------- | -------- |
| 2 × (1 × 2) | 0.1 | 1.2 × 4.3 | 2.3 × 2.2 |
| 3 × (1 × 2) | 0.1 | 2.3 × 4.3 | 3.4 × 2.2 |
| 4 × (1 × 1) | 0.5 | 3.5 × 3.5 | 3.5 × 3.5 (unchanged) |

## Acceptance Criteria
- [ ] With no `opts.target_width`, the wrap decision is made on squareness of the
      running atlas extent; the `sqrt(Σ (w+pad)(h+pad))` estimate is removed from
      the default path, so no estimate/measurement padding mismatch remains.
- [ ] Two 1 × 2 charts with `padding = 0.1` pack side by side to a 2.3 × 2.2
      extent (regression test for the reported case), and
      `examples/src/MultiChartFlatten.cpp` reports a near-square normalized
      extent instead of `[0, 0] -> [0.27907, 1]`.
- [ ] A test asserts the **default** heuristic's atlas aspect ratio on chart sets
      where near-square is achievable (e.g. 4, 6, and 9 unit charts, with and
      without padding), closing the gap noted in #104 that the existing
      pinned-extent tests set `target_width` and therefore do not constrain the
      heuristic.
- [ ] Adding padding does not change how many charts land on a shelf for a set of
      identical charts (relative to `padding = 0`) — the specific symptom of the
      mismatch.
- [ ] An explicit `opts.target_width` still wraps exactly as documented today
      (both the no-wrap and wrap cases are pinned by tests), and a chart wider
      than the target is still placed alone at the start of a shelf rather than
      wrapping an empty shelf.
- [ ] All existing `ChartPacking` tests pass unmodified, including
      `PaddingSeparatesChartsInSingleRow` and `PaddingSurroundsChartsAtPerimeter`;
      `PaddingSurroundsChartsAtPerimeter` is confirmed to still span ≥ 2 shelves
      under the new decision (its comment states it exercises multiple shelves) —
      if the new layout puts its 5 charts in one row, the chart set is adjusted so
      the multi-shelf perimeter margins remain covered.
- [ ] `PackOptions::target_width` Doxygen no longer documents a `sqrt(total area)`
      default and instead describes the squareness-based wrap and the override's
      exact (unpadded right edge vs. `target_width`) semantics; `padding` docs
      remain accurate.
- [ ] `PackCharts`'s documented `O(n log n)` complexity is unchanged (the new test
      is O(1) per chart).
- [ ] `git clang-format` clean; single header regenerated via
      `thirdparty/amalgamate/amalgamate.py`.

## Out of Scope
- Changing the semantics of an explicit `opts.target_width` (it keeps measuring
  the unpadded right edge against the given width; only its documentation is
  clarified).
- Direction 1 from #104 (repairing the sqrt-area estimate's gutter budget) — the
  squareness test removes the need for a width estimate at all.
- Replacing shelf packing with skyline/MaxRects-style placement, or changing the
  descending-height chart ordering.
- Packing-efficiency or fill-ratio reporting in the examples or benchmarks.
- Per-chart scaling; `normalize` remains a single global uniform scale.

## Dependencies
- **F2 (#18, PR #99) must merge first** — `include/OpenABF/ChartPacking.hpp` and
  its tests land there. Work starts from `develop` after that merge, not from the
  `f2-multi-chart-packing` branch.
- B10 (#103) also rides on PR #99 but is otherwise independent of P7.
