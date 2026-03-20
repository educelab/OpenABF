# P2 Implementation Plan

## Phase 1 — Fix (ABF)
Replace in `ABF::Compute`:
```cpp
// Before
std::map<std::size_t, std::size_t> vIdx2vIntIdx;
std::size_t newIdx{0};
for (const auto& v : mesh->vertices_interior()) {
    vIdx2vIntIdx[v->idx] = newIdx++;
}
// lookup: vIdx2vIntIdx.at(v->idx)

// After
std::vector<std::size_t> vIdx2vIntIdx(mesh->num_vertices(), 0);
std::size_t newIdx{0};
for (const auto& v : mesh->vertices_interior()) {
    vIdx2vIntIdx[v->idx] = newIdx++;
}
// lookup: vIdx2vIntIdx[v->idx]
```

## Phase 2 — Fix (ABFPlusPlus)
Apply the same change to `ABFPlusPlus::Compute`.

## Phase 3 — Verify
- Run `ctest`
