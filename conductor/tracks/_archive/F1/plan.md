# F1 Implementation Plan

See [design.md](design.md) for algorithm details, API design, and rationale.

## Phase 1 — Test Fixtures

Add mesh builders to `tests/src/Utils.hpp`.

1. `ConstructHemisphere<MeshType>(rings, sectors)`: generate a UV hemisphere
   (half-sphere, open at the equator). Vertex (0,0) is the pole; the last ring
   is the equatorial boundary. `rings` latitude bands × `sectors` longitude
   slices → `rings * sectors + 1` vertices, `2 * (rings-1) * sectors + sectors`
   faces. Default suggestion: `ConstructHemisphere(12, 24)` → 289 vertices,
   552 faces. Non-zero Gaussian curvature everywhere; the conformal map is
   stereographic projection, giving an analytically verifiable reference.

2. `ConstructWavySurface<MeshType>(rows, cols)`: a grid whose z-coordinates are
   displaced by `z = 0.3 * sin(2πx / (cols-1)) * cos(2πy / (rows-1))`. Same
   connectivity as `ConstructGrid` but with spatially varying curvature and a
   parameterization that is *not* already conformal. Default suggestion:
   `ConstructWavySurface(20, 20)` → 400 vertices, 722 faces.

## Phase 2 — Tests

Write failing tests in `tests/src/TestParameterization.cpp` for the public API.

1. `HLSCM_Pyramid`: on the 4-vertex pyramid, HLSCM produces the same UV output
   as `AngleBasedLSCM` (mesh too small for hierarchy → single-level fallback)
2. `HLSCM_ExplicitPins`: `HierarchicalLSCM::Compute(mesh, 0, 1)` matches the
   auto-pin result when vertex 0 and 1 are the default pins
3. `HLSCM_SetPinnedVertices`: instance API `setPinnedVertices()` + `compute()`
   matches the static `Compute(mesh, pin0, pin1)` overload
4. `HLSCM_Double`: double-precision test on the pyramid (matches
   `AngleBasedLSCM_Double` expectations)
5. `HLSCM_ABFPlusPlus`: `ABFPlusPlus::Compute` then `HierarchicalLSCM::Compute`
   produces a valid parameterization on the pyramid
6. `HLSCM_Hemisphere`: on a 12×24 hemisphere (~289 vertices), HLSCM and
   `AngleBasedLSCM` produce UV coordinates within tolerance (`1e-4` for float).
   Also verify all UVs finite, z = 0, no triangle flips.
7. `HLSCM_WavySurface`: on a 20×20 wavy surface (400 vertices), HLSCM and
   `AngleBasedLSCM` produce UV coordinates within tolerance. Exercises
   multi-level hierarchy with non-trivial curvature variation.

## Phase 3 — Internal Decimation Infrastructure

Build a lightweight mesh decimation engine in `detail::hlscm` namespace inside
`HierarchicalLSCM.hpp`. No modifications to `HalfEdgeMesh`.

1. Flat adjacency structure: copy vertex positions and face connectivity from
   `HalfEdgeMesh` into `std::vector`s; build vertex→face and edge→face maps
2. QEM quadric: compute initial Garland-Heckbert quadric per vertex from
   incident face planes; implement `Quadric<T>` with `evaluate()` and `+=`
3. Half-edge collapse: given edge `(v_remove, v_keep)`, remove incident faces,
   redirect remaining faces, merge quadrics, record collapse for prolongation
4. Priority-queue-driven greedy decimation with validity checks:
   - Link condition (no non-manifold topology post-collapse)
   - Boundary vertex collapse allowed, with guard against two boundary
     vertices collapsing via an interior edge (would create non-manifold)
   - Minimum angle threshold (10°) on all post-collapse faces
   - Normal-flip and degenerate-face rejection
   - Pinned vertices uncollapsible
5. Hierarchy builder: given a target level ratio, repeatedly collapse cheapest
   valid edge; at each level boundary (vertex count crosses a ratio threshold),
   snapshot the surviving vertices/faces and the collapse records
6. Level mesh construction: build a `HalfEdgeMesh<T>` per level from the
   snapshot using `insert_vertex()` / `insert_faces()`, maintaining an index
   map from fine vertices to coarse vertices

## Phase 4 — HLSCM Core Solver

1. Pin selection: reuse `AngleBasedLSCM`'s boundary-walk logic for auto pin;
   mark pin vertices as uncollapsible before decimation
2. Build hierarchy (Phase 2 output): sequence of `HalfEdgeMesh<T>` levels +
   per-level vertex maps + collapse records with barycentric coordinates
3. Coarsest-level solve: compute face angles via `ComputeMeshAngles()`, build
   LSCM system (A, b) using the Lévy et al. Eq. 10 formulation, solve with
   `Solver::solve(b_dense)` (no initial guess)
4. Prolongation operator: for each level k → k+1, assign surviving-vertex UVs
   directly; assign new-vertex UVs via barycentric interpolation from the
   collapse record's containing triangle
5. Finer-level solve: build LSCM system at level k+1 (with mesh angles from
   `ComputeMeshAngles()` for intermediate levels, or the input mesh's existing
   angles at the finest level), extract free-vertex UVs from prolongated
   solution as initial guess x₀, solve with `Solver::solveWithGuess(b_dense, x0)`
6. Transfer final UVs from the finest internal mesh back to the input mesh's
   vertex positions

## Phase 5 — API & Integration

1. Static `Compute(mesh)` with auto pin selection
2. Static `Compute(mesh, pin0, pin1)` with explicit pins
3. Instance `compute(mesh)` respecting `pinnedVertices_`, `levelRatio_`,
   `minCoarseVertices_`
4. `#include "OpenABF/HierarchicalLSCM.hpp"` in `OpenABF.hpp` (after
   `AngleBasedLSCM.hpp`)
5. Regenerate amalgamated header:
   `python3 thirdparty/amalgamate/amalgamate.py -c single_include.json -s .`

## Phase 6 — Verify

- `ctest` — all existing + new tests pass
- Examples compile
- `git clang-format`
- HLSCM on pyramid matches LSCM output exactly (single-level equivalence)
- HLSCM on hemisphere and wavy surface matches LSCM within tolerance
- Timing comparison: HLSCM vs AngleBasedLSCM on `ConstructWavySurface(50, 50)`
  logged to stdout (informational, not gated)
