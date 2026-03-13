# P3 Implementation Plan

## Phase 1 — Fix
In `InitializeAnglesAndWeights`, compute angle_sum and update phi in a single pass:
```cpp
for (auto& v : m->vertices_interior()) {
    auto wheel = v->wheel();
    T angle_sum = T(0);
    for (const auto& e : wheel) { angle_sum += e->beta; }
    for (auto& e : wheel) {
        e->phi *= 2 * PI<T> / angle_sum;
        e->weight = T(1) / (e->phi * e->phi);
    }
}
```

## Phase 2 — Verify
- Run `ctest`
