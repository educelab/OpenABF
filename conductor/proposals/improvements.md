# OpenABF Improvement Proposals

Findings from a code review across all headers, tests, and the GitHub issue tracker.
Issues are categorized by severity and type.

---

## Bugs

### B1 — ABF lambda update index is wrong for meshes with 2+ interior vertices

**File:** `include/OpenABF/ABF.hpp` (~line 398)

The interior-vertex lambda update loop uses `idx + intIdx` where both `idx` and `intIdx`
increment each iteration, producing incorrect offsets into the delta vector:

```cpp
for (auto& v : mesh->vertices_interior()) {
    auto intIdx = vIdx2vIntIdx.at(v->idx);
    v->lambda_plan += delta(idx + intIdx, 0);       // wrong for i > 0
    v->lambda_len  += delta(idx + vIntCnt + intIdx, 0);
    idx++;
}
```

For interior vertex `i`, at loop entry `idx = edgeCnt + faceCnt + i` and `intIdx = i`,
so the offset is `edgeCnt + faceCnt + 2*i` instead of the correct
`edgeCnt + faceCnt + i`. The corresponding ABFPlusPlus code gets this right using
`faceCnt + intIdx` directly. The pyramid used in the unit test has only one interior
vertex (`vIntCnt = 1`), which masks the bug since `i` and `intIdx` are both 0.

**Fix:** Replace `idx + intIdx` with a fixed base: `auto base = edgeCnt + faceCnt` and
use `base + intIdx` for lambda_plan, `base + vIntCnt + intIdx` for lambda_len.

---

### B2 — `ABFPlusPlus` LambdaStarInv inversion uses hardcoded `1.F` (float) regardless of scalar type `T`

**File:** `include/OpenABF/ABFPlusPlus.hpp` (line 187)

```cpp
it.valueRef() = 1.F / it.value();
```

When `T = double`, this truncates the result to single precision. Should be `T(1) / it.value()`.

---

### B3 — `Vec` reverse iterators have wrong return type

**File:** `include/OpenABF/Vec.hpp` (lines 110–127)

`rbegin()`, `rend()`, `crbegin()`, `crend()` are declared to return `iterator` /
`const_iterator` but should return `reverse_iterator` / `const_reverse_iterator`.
This causes silent type mismatches when used with reverse iteration.

---

### B4 — `Vec::operator*` and `operator/` are inconsistent with `operator*=` and `operator/=`

**File:** `include/OpenABF/Vec.hpp`

`operator*=(T2)` and `operator/=(T2)` require an arithmetic scalar, but `operator*(Vec, Vector)`
calls `*=` with a general `Vector`, which will fail to compile or produce incorrect behavior
when a non-scalar vector is passed. The binary `*` and `/` operators should either be
constrained to scalars (like `*=`) or the `*=` should be generalized to element-wise
operations. Currently the operators are misleading.

---

### B5 — `is_file_type` in `MeshIOFormats.hpp` accesses `ext[0]` without bounds check

**File:** `include/OpenABF/MeshIOFormats.hpp` (line 22)

```cpp
auto ext = path.extension().string();
if (ext[0] == '.') {  // UB if ext is empty (no extension)
```

If a path has no extension, `ext` is an empty string and `ext[0]` is undefined behavior.

---

### B6 — PLY reader: `vmap` array may be uninitialized if x/y/z properties are absent

**File:** `include/OpenABF/MeshIOFormats.hpp`

`std::array<std::size_t, 3> vmap{}` is value-initialized to `{0, 0, 0}`. If the PLY header
is missing one or more of `x`, `y`, `z` properties, those indices stay at 0 and the
reader silently reads the wrong column. Should validate that all three are found.

---

### B7 — `PlanGrad` allocates and discards the wheel vector

**File:** `include/OpenABF/ABF.hpp` (line 98)

```cpp
auto edges = v->wheel();  // allocated but never used
T g = -2 * PI<T>;
for (const auto& e : v->wheel()) {  // second allocation
```

The first `v->wheel()` result is unused. Minor bug / dead code, but wastes a heap
allocation on every call.

---

## Performance

### P1 — All mesh iteration methods return `std::vector` by value

**File:** `include/OpenABF/HalfEdgeMesh.hpp`
**Also:** GitHub issue #4

`vertices()`, `edges()`, `faces()`, `vertices_interior()`, `vertices_boundary()`,
`wheel()`, `boundaries()`, `outgoing_edges()`, `connected_components()` all allocate
and return `std::vector` by value. The ABF solver calls these in tight iteration loops
(per-iteration of the outer gradient loop), causing O(n) heap allocations per solver
step that are immediately discarded.

Consider replacing with lazy range views (C++20 `std::ranges`, or a lightweight custom
`IterableView`) or cached pre-computed lists invalidated on mutation.

---

### P2 — `std::map` used for interior-vertex index lookups in hot paths

**File:** `include/OpenABF/ABF.hpp`, `include/OpenABF/ABFPlusPlus.hpp`

`vIdx2vIntIdx` is a `std::map<std::size_t, std::size_t>` rebuilt every solver iteration
with O(log n) lookup. Since vertex indices are dense (0..N-1), a flat `std::vector`
with direct indexing would be O(1) and avoid per-lookup tree traversal.

---

### P3 — `InitializeAnglesAndWeights` calls `v->wheel()` redundantly

**File:** `include/OpenABF/ABF.hpp` (line 73)

```cpp
auto wheel = v->wheel();  // first allocation
auto angle_sum = std::accumulate(wheel.begin(), wheel.end(), T(0), ...);
for (auto& e : wheel) { ... }  // one pass could do both
```

This is correct but accumulate and update could be two passes over the same vector.
With P1 fixed, this would also benefit automatically.

---

## API / Usability

### A1 — Convergence tolerance is hardcoded

**Files:** `include/OpenABF/ABF.hpp`, `include/OpenABF/ABFPlusPlus.hpp`

The gradient threshold `0.001` and gradient-delta threshold `0.001` are hardcoded in
the solver loop. Users working with very large or very small meshes, or needing different
precision/speed trade-offs, cannot tune these. A `setGradientThreshold(T)` method and/or
constructor parameter should be added alongside the existing `setMaxIterations`.

---

### A2 — `FacePtr` iteration requires awkward `*face` dereference

**File:** `include/OpenABF/HalfEdgeMesh.hpp`
**Also:** GitHub issue #4

Range-based iteration over a face pointer requires `for (const auto& e : *face)`.
The proposed fix from the issue tracker is a `face->edges()` method that returns an
iterable. A named method makes intent clearer and avoids the dereference of a
`shared_ptr` in user code.

---

### A3 — `insert_face` (single-call form) boundary update behavior is undocumented and inconsistent

**File:** `include/OpenABF/HalfEdgeMesh.hpp`

The variadic `insert_face(Args... args)` calls `insert_face_(...)` but does NOT call
`update_boundary()`, matching the behavior of the vector-based `insert_face(Vector&&)`.
However, users expecting a single face insertion to be a complete operation may be
surprised. The documentation on the vector form says "does not update boundary" but
the variadic form omits this caveat. Consistency and documentation should be improved.

---

### A4 — `gradient()` is not `[[nodiscard]]` but `iterations()` is

**Files:** `include/OpenABF/ABF.hpp`, `include/OpenABF/ABFPlusPlus.hpp`

Both methods return state that is only meaningful after calling `compute()`. Applying
`[[nodiscard]]` consistently to both, and ideally documenting their initial values,
would help users avoid accidentally reading uninitialized state.

---

### A5 — No configurable pinned edge selection for LSCM

**File:** `include/OpenABF/AngleBasedLSCM.hpp`
**Also:** GitHub issue #2

The pinned edge is always the first boundary edge found. Blender's LSCM uses a
heuristic to select an edge that balances the chart. An optional callback or strategy
object for pinned-edge selection would allow callers to use better heuristics without
forking the implementation.

---

### A6 — AB matrix construction in AngleBasedLSCM is hard to follow

**File:** `include/OpenABF/AngleBasedLSCM.hpp`
**Also:** GitHub issue #3

The A and B matrix assembly in `AngleBasedLSCM::Compute` is a large block of
repetitive per-vertex conditionals. Extracting a helper lambda or function for the
per-vertex contribution would improve readability and testability.

---

## Missing Features

### F1 — HLSCM (Hierarchical LSCM)

**Also:** GitHub issue #5 (references the Levy et al. paper)

No hierarchical LSCM implementation exists. The referenced paper describes an approach
that produces better parameterizations on high-curvature regions by hierarchically
refining the parameterization. This would be a significant addition.

---

### F2 — Multi-chart packing

A key project goal is "multi-chart flattening and packing." There is no chart packing
API. After splitting a mesh with `split_path` into multiple connected components and
running LSCM per-component, there is no utility to pack charts into a [0,1]² UV
atlas. A basic rectangle-packing or shelf-packing algorithm for UV charts would
address this gap.

---

### F3 — Multi-component extraction and parameterization pipeline

There is no way to feed a mesh with multiple connected components (e.g. after one or
more `split_path` calls) directly into the parameterization pipeline. ABF and LSCM
both assume a single connected component.

**Proposed design — two layers:**

**Low-level extractor** (free function, consistent with `FindEdgePath`, `ComputeMeshAngles`):

```cpp
template <class MeshType>
auto ExtractConnectedComponents(const typename MeshType::Pointer& mesh)
    -> std::vector<std::pair<typename MeshType::Pointer,
                             std::vector<std::size_t>>>;
// pair.second[extracted_idx] == original_idx
```

Each extracted mesh is a full deep copy of the CC's vertices, edges, faces, and all
traits (including pre-computed `alpha`/`beta`/`phi` values if `ComputeMeshAngles` was
run before extraction). Vertex indices are re-densified (0..N-1) in the new mesh, and
the back-map vector provides the correspondence to the original mesh for writing UV
coordinates back or building a texture atlas.

The extracted meshes are fully independent, so callers can parallelize with any
threading model (`std::async`, OpenMP, TBB) without the library mandating a choice.

**High-level convenience helper:**

```cpp
template <class AngleOptimizer, class Parameterizer, class MeshType>
void ParameterizeConnectedComponents(typename MeshType::Pointer& mesh);
```

Calls `ExtractConnectedComponents`, runs the `AngleOptimizer` + `Parameterizer`
pipeline on each CC in sequence, and writes the UV results back to the original mesh's
vertex positions. This covers the common case where the caller doesn't need to
inspect intermediate per-CC results.

**Typical end-to-end workflow this enables:**

```
split_path → ExtractConnectedComponents → [ABF++ + LSCM] per CC → pack into atlas
```

**Implementation notes:**
- Trait copying works naturally through the existing copy constructors
  (`Edge(const Edge& rhs) : EdgeTraits(rhs) {}`, etc.).
- The extractor can build on the existing `connected_components()` method internally.
- Should be a new header (`HalfEdgeMeshUtils.hpp` or similar) rather than added to
  the already large `HalfEdgeMesh.hpp`.

---

### F4 — No `double`-precision test coverage

All parameterization tests use `float`. Given the numerical accuracy goal, tests
verifying that `double` produces more accurate results (e.g. tighter error bounds
against analytical solutions for known-geometry meshes) would strengthen confidence
in the library.

---

## Test Coverage Gaps

### T1 — ABF only tested on a 4-vertex pyramid (masks bug B1)

The pyramid has exactly one interior vertex. Any bug that only manifests with 2+
interior vertices (like B1) is invisible. Tests should include:
- A mesh with ≥ 3 interior vertices (e.g. a flat grid)
- A mesh with 0 interior vertices (degenerate: should not crash)
- A mesh where ABF/ABFPlusPlus hits `maxIters` without converging

---

### T2 — No test for `FindEdgePath` on disconnected mesh

`FindEdgePath` returns an empty vector if no path exists, but this case is untested.

---

### T3 — No test for PLY/OBJ round-trip accuracy

Reading a mesh, writing it, and reading it back should produce an identical mesh.
No such round-trip test exists.

---

### T4 — No test for `Vec` reverse iterators (masks bug B3)

The incorrect return types on `rbegin()`/`rend()` are not caught because there are
no tests that use reverse iteration on `Vec`.

---

## Minor / Cleanup

### M1 — `Face::barycenter()` hardcodes `Vec<T, 3>` instead of using `Dim`

**File:** `include/OpenABF/HalfEdgeMesh.hpp`

`auto barycenter() const -> Vec<T, 3>` ignores the template parameter `Dim`.
Should be `Vec<T, Dim>`.

---

### M2 — `Edge::magnitude()` is non-const

**File:** `include/OpenABF/HalfEdgeMesh.hpp`

`auto magnitude() -> T` is non-const even though it only reads vertex positions.
This prevents calling it on const-qualified edge references and is used in `Face::area()`
(which is const), requiring a non-const path through shared_ptr.

---

### M3 — `detail::erase_if` takes container by value (shadowed by C++20 std::erase_if)

**File:** `include/OpenABF/HalfEdgeMesh.hpp`

The internal `detail::erase_if` copies the container, modifies the copy, and returns it.
Callers that don't capture the return value silently do nothing. This function is not
used anywhere in the project (grep finds no call sites). It can be removed.

---

### M4 — `operator<<` for `Vec` is defined outside the `OpenABF` namespace

**File:** `include/OpenABF/Vec.hpp` (line 274)

The streaming operator is defined after the namespace closing brace. While technically
valid (ADL will find it), it is cleaner to define it inside the namespace or as an
inline friend inside the class.

---

## Summary Table

| ID | Category | Severity | Effort | GitHub Issue |
|----|----------|----------|--------|--------------|
| B1 | Bug      | High     | Low    | —            |
| B2 | Bug      | Medium   | Trivial| —            |
| B3 | Bug      | Medium   | Low    | —            |
| B4 | Bug      | Medium   | Low    | —            |
| B5 | Bug      | Low      | Trivial| —            |
| B6 | Bug      | Low      | Low    | —            |
| B7 | Bug      | Low      | Trivial| —            |
| P1 | Perf     | High     | High   | #4           |
| P2 | Perf     | Medium   | Low    | —            |
| P3 | Perf     | Low      | Trivial| —            |
| A1 | API      | Medium   | Low    | —            |
| A2 | API      | Medium   | Low    | #4           |
| A3 | API      | Low      | Low    | —            |
| A4 | API      | Low      | Trivial| —            |
| A5 | API      | Medium   | Medium | #2           |
| A6 | API      | Low      | Medium | #3           |
| F1 | Feature  | High     | High   | #5           |
| F2 | Feature  | High     | High   | —            |
| F3 | Feature  | High     | Medium | —            |
| F4 | Feature  | Medium   | Medium | —            |
| T1 | Testing  | High     | Low    | —            |
| T2 | Testing  | Low      | Trivial| —            |
| T3 | Testing  | Low      | Low    | —            |
| T4 | Testing  | Low      | Trivial| —            |
| M1 | Cleanup  | Low      | Trivial| —            |
| M2 | Cleanup  | Low      | Trivial| —            |
| M3 | Cleanup  | Low      | Trivial| —            |
| M4 | Cleanup  | Low      | Trivial| —            |
