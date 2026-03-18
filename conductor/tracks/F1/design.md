# F1 Design — HLSCM (Hierarchical LSCM)

## Paper: Ray & Lévy, "Hierarchical Least Squares Conformal Map" (2003)

### Core Insight

Minimizers of the conformal energy are harmonic maps — barycentric coordinates
are locally preserved through the parameterization. If a vertex is inserted into
a triangle, its optimal UV equals its barycentric interpolation in that triangle.
Consequently, a conformal parameterization of a coarse mesh is already a good
approximation for a finer mesh, making cascadic multigrid highly effective.

### Algorithm

1. **Decimate** the input mesh via half-edge collapse to build a progressive mesh
   hierarchy. Level boundaries are placed at every ~10× vertex ratio.

2. **Solve** standard LSCM on the coarsest mesh (full solve).

3. **Prolongate** UV coordinates to the next finer level: each newly inserted
   vertex gets UVs via barycentric interpolation in its containing coarse
   triangle.

4. **Refine** at the finer level using conjugate gradient with the prolongated
   UVs as the initial guess. Convergence is fast because the guess is close to
   the optimum (harmonicity guarantee).

5. **Repeat** steps 3–4 until the finest level is reached.

### Decimation Details (from paper)

- Uses half-edge collapse (`hcoll`) to coarsen the mesh.
- Edge collapse priority based on geometric importance:
  - Volume difference between original and simplified mesh
  - Area difference (for boundary edges)
  - Accumulated loss of quality on merged vertices
- Validity checks: collapse must not create non-manifold topology; angular
  defect of the vertex to remove must be < 70°.
- Each `vsplit` (reverse of `hcoll`) inserts a vertex whose UV is interpolated
  barycentrically in the nearest coarse triangle.

### Performance (from paper, 2003 hardware)

| Mesh      | Triangles | LSCM (monogrid) | HLSCM (total) | Speedup |
|-----------|-----------|------------------|---------------|---------|
| Venus     | 37K       | 50s              | 18s           | 2.8×    |
| Horse     | 72K       | 390s             | 31s           | 12.6×   |
| Scanned   | 1.12M     | did not finish   | 704s          | ∞       |

---

## API Design

```cpp
template <typename T, class MeshType = HalfEdgeMesh<T>,
          class Solver = Eigen::ConjugateGradient<Eigen::SparseMatrix<T>,
                                                   Eigen::Lower | Eigen::Upper>,
          std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
class HierarchicalLSCM
{
public:
    using Mesh = MeshType;

    /** @brief Set pinned vertices (same API as AngleBasedLSCM) */
    void setPinnedVertices(std::size_t pin0Idx, std::size_t pin1Idx);

    /** @brief Vertex ratio between consecutive hierarchy levels (default: 10) */
    void setLevelRatio(std::size_t ratio);

    /** @brief Minimum vertex count at coarsest level (default: 100) */
    void setMinCoarseVertices(std::size_t count);

    /** @brief Compute with instance configuration */
    void compute(typename Mesh::Pointer& mesh) const;

    /** @brief Compute with auto pin selection (same as AngleBasedLSCM) */
    static void Compute(typename Mesh::Pointer& mesh);

    /** @brief Compute with explicit pin indices */
    static void Compute(typename Mesh::Pointer& mesh,
                        std::size_t pin0Idx, std::size_t pin1Idx);

private:
    std::optional<std::pair<std::size_t, std::size_t>> pinnedVertices_;
    std::size_t levelRatio_{10};
    std::size_t minCoarseVertices_{100};
};
```

### Key API Notes

- **Solver default**: `ConjugateGradient<SparseMatrix<T>, Lower|Upper>` which
  solves the normal equations (AᵀA x = Aᵀb) and enables OpenMP-parallelized
  SpMV. `LeastSquaresConjugateGradient` is a valid alternative but lacks the
  OpenMP benefit. Direct solvers get no benefit from the hierarchy since they
  solve from scratch every time.

- **Drop-in replacement**: `HierarchicalLSCM<T>::Compute(mesh)` can replace
  `AngleBasedLSCM<T>::Compute(mesh)` with no other changes.

- **Graceful degradation**: Meshes with fewer vertices than `minCoarseVertices`
  are solved in a single LSCM pass with no hierarchy overhead.

- **ABF compatibility**: Works with `ABF::Mesh` and `ABFPlusPlus::Mesh` as
  `MeshType`, just like `AngleBasedLSCM`.

---

## Internal Architecture

### 1. Decimation (internal, no HalfEdgeMesh modifications)

A lightweight flat-array decimation engine lives entirely within HLSCM's
implementation `detail` namespace. This avoids modifying the core `HalfEdgeMesh`
data structure.

**Data structures:**
- Vertex positions, face connectivity, vertex→face adjacency in flat
  `std::vector`s.
- Garland-Heckbert QEM quadric per vertex. Half-edge collapse cost =
  `(Q_remove + Q_keep).evaluate(v_keep)`.
- `std::priority_queue` for greedy cheapest-edge collapse.

**Validity checks per collapse:**
- No non-manifold topology (link condition; reject two boundary vertices
  collapsing via an interior edge)
- Minimum angle threshold (10°) on all post-collapse faces
- Normal-flip rejection on modified faces
- Degenerate face rejection (duplicate vertices after substitution)
- Pinned vertices are never collapsed

**Collapse records:** Each collapse stores `{v_removed, v_kept, containing_face,
barycentric_coords}` for use during prolongation.

**Level snapshots:** At each level boundary (every `levelRatio`× vertices), a new
`HalfEdgeMesh<T>` is built from surviving vertices and faces using the standard
`insert_vertex()` / `insert_faces()` construction API.

### 2. Prolongation

After solving level k, UV coordinates for level k+1:
- **Surviving vertices**: UVs inherited directly (same vertex in both levels).
- **New vertices**: barycentric interpolation in the coarse triangle recorded
  during decimation.

### 3. LSCM System Building & Solve

Each level's LSCM system uses the same formulation as
`AngleBasedLSCM::ComputeImpl` (Lévy et al. 2002, Eq. 10):

1. Select pins (same pin selection as `AngleBasedLSCM`).
2. Build the free-vertex index table.
3. Assemble sparse matrices A (free) and bFree (pinned) from edge angles.
4. Compute RHS: `b = -bFree * bFixed`.
5. Solve:
   - **Coarsest level**: `solver.solve(b_dense)` (no initial guess).
   - **Finer levels**: `solver.solveWithGuess(b_dense, x0)` where x0 is the
     prolongated UV vector.

Note: `solveWithGuess` requires dense matrices, so sparse b is converted to
dense before solving.

**Edge angles**: Coarse levels compute angles from 3D geometry via
`ComputeMeshAngles()`. The finest level uses the input mesh's existing angles
(which may be ABF-optimized).

### 4. Pin Mapping Across Levels

The two pinned vertices must exist at every level. During decimation, pinned
vertices are marked as uncollapsible. The pin indices are mapped from the fine
mesh to each coarse mesh using the vertex mapping.

---

## Design Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Decimation location | Internal to HLSCM | Avoids modifying HalfEdgeMesh; self-contained |
| Error metric | QEM (Garland-Heckbert) | Standard, good geometric fidelity |
| Level ratio | 10× (configurable) | Paper recommendation |
| Default solver | CG with `Lower\|Upper` | OpenMP-parallelized SpMV; best multi-thread perf |
| CG tolerance | 1e-8 | Sufficient for UV; enables effective warm-start |
| Coarse mesh type | `HalfEdgeMesh<T>` | Simpler than user's MeshType; internal only |
| ABF interaction | Geometry angles at coarse, mesh angles at finest | Preserves ABF optimization |
| Boundary handling | Allow collapse with non-manifold guards | Matches paper; enables deep hierarchy on open meshes |
| Pin handling | Pins marked uncollapsible | Ensures pins exist at all levels |

## Future Considerations

- **A7 (multi-pin)**: Once multi-pin support lands in `AngleBasedLSCM`, HLSCM
  can forward a `PinMap` to the system builder. All pins would be marked
  uncollapsible during decimation.
- **LOD output**: The hierarchy could be exposed for progressive mesh rendering
  (paper §4.1), but this is beyond F1 scope.
- **Extracting decimation**: If F2/F3 need mesh simplification, the internal
  decimation engine could be extracted to a shared utility.
