#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <optional>
#include <queue>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCore>

#include "OpenABF/AngleBasedLSCM.hpp"
#include "OpenABF/Exceptions.hpp"
#include "OpenABF/HalfEdgeMesh.hpp"
#include "OpenABF/Math.hpp"

namespace OpenABF
{

namespace detail
{
/**
 * @brief Implementation details for HierarchicalLSCM
 */
namespace hlscm
{

/**
 * @brief Symmetric 4×4 quadric matrix for QEM error metric (Garland-Heckbert)
 */
template <typename T>
struct Quadric {
    /** Upper triangle stored row-major: a00 a01 a02 a03 a11 a12 a13 a22 a23 a33 */
    std::array<T, 10> q{};

    Quadric() = default;

    /** Construct from plane equation ax + by + cz + d = 0 */
    Quadric(T a, T b, T c, T d)
        : q{a * a, a * b, a * c, a * d, b * b, b * c, b * d, c * c, c * d, d * d}
    {
    }

    /** In-place accumulation of another quadric */
    auto operator+=(const Quadric& o) -> Quadric&
    {
        for (std::size_t i = 0; i < 10; ++i) {
            q[i] += o.q[i];
        }
        return *this;
    }

    /** Quadric addition */
    friend auto operator+(Quadric a, const Quadric& b) -> Quadric { return a += b; }

    /** Evaluate quadric error at point (x, y, z) */
    auto evaluate(T x, T y, T z) const -> T
    {
        // v^T Q v where Q is the symmetric 4x4 matrix, v = (x, y, z, 1)
        return q[0] * x * x + T(2) * q[1] * x * y + T(2) * q[2] * x * z + T(2) * q[3] * x +
               q[4] * y * y + T(2) * q[5] * y * z + T(2) * q[6] * y + q[7] * z * z +
               T(2) * q[8] * z + q[9];
    }
};

/**
 * @brief Record of a single half-edge collapse for prolongation
 */
template <typename T>
struct CollapseRecord {
    /** Index of the removed vertex (in the original/fine mesh) */
    std::size_t v_removed;
    /** Index of the kept vertex (in the original/fine mesh) */
    std::size_t v_kept;
    /** Post-collapse triangle containing v_removed (original vertex indices) */
    std::array<std::size_t, 3> containing_tri;
    /** Barycentric coordinates of v_removed in containing_tri */
    std::array<T, 3> bary;
};

/**
 * @brief A level in the mesh hierarchy
 */
template <typename T>
struct HierarchyLevel {
    /** Vertex positions (indexed by level-local index) */
    std::vector<Vec<T, 3>> positions;
    /** Face connectivity (each face is 3 level-local indices) */
    std::vector<std::array<std::size_t, 3>> faces;
    /** Map from level-local vertex index to original (finest) vertex index */
    std::vector<std::size_t> local_to_original;
    /** Map from original vertex index to level-local index. Sized to the
     *  finest-mesh vertex count; vertices not present at this level hold
     *  `std::nullopt`.
     */
    std::vector<std::optional<std::size_t>> original_to_local;
};

/**
 * @brief Lightweight flat-array mesh for decimation
 *
 * Copies vertex positions and face connectivity from a HalfEdgeMesh into
 * flat vectors, builds adjacency structures, and supports half-edge collapse.
 */
template <typename T>
class DecimationMesh
{
public:
    /** Build from a HalfEdgeMesh */
    template <class MeshPtr>
    void build(const MeshPtr& mesh, const std::vector<std::size_t>& pinIndices)
    {
        auto nv = mesh->num_vertices();
        auto nf = mesh->num_faces();

        positions_.resize(nv);
        alive_.assign(nv, true);
        is_boundary_.assign(nv, false);
        is_pinned_.assign(nv, false);
        quadrics_.resize(nv);

        for (auto idx : pinIndices) {
            is_pinned_[idx] = true;
        }

        for (const auto& v : mesh->vertices()) {
            positions_[v->idx] = v->pos;
            is_boundary_[v->idx] = v->is_boundary();
        }

        faces_.reserve(nf);
        face_alive_.reserve(nf);
        vert_faces_.resize(nv);

        for (const auto& f : mesh->faces()) {
            auto e0 = f->head;
            auto e1 = e0->next;
            auto e2 = e1->next;
            std::array<std::size_t, 3> tri{e0->vertex->idx, e1->vertex->idx, e2->vertex->idx};
            auto fi = faces_.size();
            faces_.push_back(tri);
            face_alive_.push_back(true);
            vert_faces_[tri[0]].push_back(fi);
            vert_faces_[tri[1]].push_back(fi);
            vert_faces_[tri[2]].push_back(fi);
        }

        num_alive_verts_ = nv;
        num_alive_faces_ = nf;

        compute_quadrics_();
        build_edges_();
    }

    /** Get number of alive vertices */
    [[nodiscard]] auto num_alive_verts() const -> std::size_t { return num_alive_verts_; }

    /** Get number of alive faces */
    [[nodiscard]] auto num_alive_faces() const -> std::size_t { return num_alive_faces_; }

    /**
     * @brief Try to collapse edge (vRemove → vKeep), returning a collapse record
     *
     * Returns nullopt if the collapse is invalid.
     *
     * @param outKeepNbrs  Optional out-vector. On a successful collapse it is
     *                     overwritten with vKeep's post-collapse neighbor
     *                     vertex indices (sorted, unique). Passing a caller-
     *                     owned vector here lets the PQ-update path reuse the
     *                     allocation across collapses instead of letting
     *                     `vertex_neighbors()` allocate a fresh vector each
     *                     time.
     */
    auto try_collapse(std::size_t vRemove, std::size_t vKeep,
                      std::vector<std::size_t>* outKeepNbrs = nullptr)
        -> std::optional<CollapseRecord<T>>
    {
        if (!alive_[vRemove] || !alive_[vKeep]) {
            return std::nullopt;
        }
        if (is_pinned_[vRemove]) {
            return std::nullopt;
        }

        // Find shared faces (will be removed) and vRemove-only faces (will be updated)
        // Reuse scratch storage
        auto& sharedFaces = scratch_shared_;
        auto& removeFaces = scratch_remove_;
        sharedFaces.clear();
        removeFaces.clear();

        for (auto fi : vert_faces_[vRemove]) {
            if (!face_alive_[fi])
                continue;
            bool hasKeep = false;
            for (auto vi : faces_[fi]) {
                if (vi == vKeep) {
                    hasKeep = true;
                    break;
                }
            }
            if (hasKeep)
                sharedFaces.push_back(fi);
            else
                removeFaces.push_back(fi);
        }

        // Interior edges have 2 shared faces; boundary edges have 1.
        if (sharedFaces.empty() || sharedFaces.size() > 2)
            return std::nullopt;

        // Reject collapse of two boundary vertices via an interior edge:
        // vKeep would inherit two disconnected boundary fans → non-manifold.
        if (is_boundary_[vRemove] && is_boundary_[vKeep] && sharedFaces.size() == 2) {
            return std::nullopt;
        }

        // Link condition: collect sorted unique neighbors of vRemove and vKeep
        auto fillSortedNeighbors = [&](std::size_t v, std::vector<std::size_t>& out) {
            out.clear();
            for (auto fi : vert_faces_[v]) {
                if (!face_alive_[fi])
                    continue;
                for (auto vi : faces_[fi]) {
                    if (vi != v)
                        out.push_back(vi);
                }
            }
            std::sort(out.begin(), out.end());
            out.erase(std::unique(out.begin(), out.end()), out.end());
        };

        fillSortedNeighbors(vRemove, scratch_nbrs_a_);
        fillSortedNeighbors(vKeep, scratch_nbrs_b_);

        // Shared neighbors (excluding vKeep/vRemove from each other's sets)
        scratch_shared_nbrs_.clear();
        for (auto v : scratch_nbrs_a_) {
            if (v != vKeep &&
                std::binary_search(scratch_nbrs_b_.begin(), scratch_nbrs_b_.end(), v)) {
                scratch_shared_nbrs_.push_back(v);
            }
        }
        // scratch_shared_nbrs_ is already sorted since scratch_nbrs_a_ is sorted

        // Expected shared: opposite vertices of the shared faces (at most 2 entries)
        std::array<std::size_t, 2> expectedShared{};
        std::size_t numExpected = 0;
        for (auto fi : sharedFaces) {
            for (auto vi : faces_[fi]) {
                if (vi != vRemove && vi != vKeep)
                    expectedShared[numExpected++] = vi;
            }
        }
        if (numExpected > 1 && expectedShared[0] > expectedShared[1])
            std::swap(expectedShared[0], expectedShared[1]);

        if (scratch_shared_nbrs_.size() != numExpected)
            return std::nullopt;
        for (std::size_t i = 0; i < numExpected; ++i) {
            if (scratch_shared_nbrs_[i] != expectedShared[i])
                return std::nullopt;
        }

        // Minimum angle threshold (radians) — reject collapses that would
        // create triangles with any angle below this.  Paper uses ~10°.
        constexpr T minAngle = PI<T> / T(18);  // 10°

        // Validate all faces that will exist around vKeep after collapse:
        // removeFaces (with vRemove→vKeep substitution) must not flip or
        // degenerate, and ALL surviving faces incident to vKeep must
        // maintain a minimum angle above the threshold.
        //
        // Collect the full set of post-collapse faces incident to vKeep.
        auto& postFaces = scratch_post_faces_;
        postFaces.clear();
        // Existing vKeep faces (excluding shared faces which will be removed)
        auto isInShared = [&](std::size_t fi) {
            return std::find(sharedFaces.begin(), sharedFaces.end(), fi) != sharedFaces.end();
        };
        for (auto fi : vert_faces_[vKeep]) {
            if (!face_alive_[fi] || isInShared(fi)) {
                continue;
            }
            postFaces.push_back(faces_[fi]);
        }
        // removeFaces with vRemove→vKeep substitution
        for (auto fi : removeFaces) {
            std::array<std::size_t, 3> newTri = faces_[fi];
            for (auto& vi : newTri) {
                if (vi == vRemove) {
                    vi = vKeep;
                }
            }
            // Reject if substitution creates a degenerate face
            if (newTri[0] == newTri[1] || newTri[1] == newTri[2] || newTri[0] == newTri[2]) {
                return std::nullopt;
            }
            postFaces.push_back(newTri);

            // Also check for normal flip on the modified faces
            auto& op0 = positions_[faces_[fi][0]];
            auto& op1 = positions_[faces_[fi][1]];
            auto& op2 = positions_[faces_[fi][2]];
            auto oldNormal = cross(op1 - op0, op2 - op0);
            auto newNormal = cross(positions_[newTri[1]] - positions_[newTri[0]],
                                   positions_[newTri[2]] - positions_[newTri[0]]);
            if (dot(oldNormal, newNormal) < T(0)) {
                return std::nullopt;
            }
        }

        // Check all post-collapse faces for minimum angle
        for (auto& tri : postFaces) {
            auto& p0 = positions_[tri[0]];
            auto& p1 = positions_[tri[1]];
            auto& p2 = positions_[tri[2]];
            auto e01 = p1 - p0;
            auto e02 = p2 - p0;
            auto e12 = p2 - p1;
            auto l01 = norm(e01);
            auto l02 = norm(e02);
            auto l12 = norm(e12);
            if (l01 == T(0) || l02 == T(0) || l12 == T(0)) {
                return std::nullopt;
            }
            // Clamp acos argument to [-1,1] for numerical safety
            auto clampedAngle = [](T cosVal) -> T {
                return std::acos(std::max(T(-1), std::min(T(1), cosVal)));
            };
            T a0 = clampedAngle(dot(e01, e02) / (l01 * l02));
            T a1 = clampedAngle(dot(p0 - p1, e12) / (l01 * l12));
            T a2 = PI<T> - a0 - a1;
            if (a0 < minAngle || a1 < minAngle || a2 < minAngle) {
                return std::nullopt;
            }
        }

        // Build collapse record with barycentric coordinates
        // After collapse, face (vRemove, vA, vB) becomes (vKeep, vA, vB).
        // Store bary coords of v_removed's position in the post-collapse triangle.
        CollapseRecord<T> record;
        record.v_removed = vRemove;
        record.v_kept = vKeep;

        if (!removeFaces.empty()) {
            // Use the first surviving face; its post-collapse vertices are
            // (vKeep, vA, vB) where vA and vB are the non-vRemove vertices.
            auto fi = removeFaces[0];
            std::array<std::size_t, 3> postTri;
            postTri[0] = vKeep;
            std::size_t slot = 1;
            for (auto vi : faces_[fi]) {
                if (vi != vRemove) {
                    postTri[slot++] = vi;
                }
            }
            record.containing_tri = postTri;
            record.bary = computeBarycentric_(positions_[postTri[0]], positions_[postTri[1]],
                                              positions_[postTri[2]], positions_[vRemove]);
        } else {
            // Edge case: all faces are shared — vertex collapses directly onto vKeep
            record.containing_tri = {vKeep, vKeep, vKeep};
            record.bary = {T(1), T(0), T(0)};
        }

        // Execute collapse
        // Kill shared faces
        for (auto fi : sharedFaces) {
            face_alive_[fi] = false;
            num_alive_faces_--;
        }

        // Update removeFaces: replace vRemove with vKeep
        for (auto fi : removeFaces) {
            for (auto& vi : faces_[fi]) {
                if (vi == vRemove) {
                    vi = vKeep;
                }
            }
            vert_faces_[vKeep].push_back(fi);
        }

        // Mark vRemove as dead
        alive_[vRemove] = false;
        num_alive_verts_--;

        // Propagate boundary status: if vRemove was on the boundary,
        // vKeep inherits it (it now sits on the mesh boundary).
        if (is_boundary_[vRemove]) {
            is_boundary_[vKeep] = true;
        }

        // Merge quadrics
        quadrics_[vKeep] += quadrics_[vRemove];

        // Compact dead face indices from vKeep's adjacency list
        auto& vkFaces = vert_faces_[vKeep];
        vkFaces.erase(std::remove_if(vkFaces.begin(), vkFaces.end(),
                                     [this](std::size_t fi) { return !face_alive_[fi]; }),
                      vkFaces.end());

        // Fill caller-owned post-collapse neighbor vector (reusing its
        // allocation across calls). Equivalent to vertex_neighbors(vKeep) but
        // without an additional heap allocation.
        if (outKeepNbrs != nullptr) {
            outKeepNbrs->clear();
            for (auto fi : vkFaces) {
                for (auto vi : faces_[fi]) {
                    if (vi != vKeep && alive_[vi]) {
                        outKeepNbrs->push_back(vi);
                    }
                }
            }
            std::sort(outKeepNbrs->begin(), outKeepNbrs->end());
            outKeepNbrs->erase(std::unique(outKeepNbrs->begin(), outKeepNbrs->end()),
                               outKeepNbrs->end());
        }

        return record;
    }

    /** Compute collapse cost for edge (v0 → v1): Q_merged evaluated at v1 */
    [[nodiscard]] auto collapse_cost(std::size_t v0, std::size_t v1) const -> T
    {
        auto Q = quadrics_[v0] + quadrics_[v1];
        auto& p = positions_[v1];
        return Q.evaluate(p[0], p[1], p[2]);
    }

    /** Check if a vertex is alive */
    [[nodiscard]] auto is_alive(std::size_t v) const -> bool { return alive_[v]; }

    /** Check if a vertex is collapsible (not pinned, alive) */
    [[nodiscard]] auto is_collapsible(std::size_t v) const -> bool
    {
        return alive_[v] && !is_pinned_[v];
    }

    /** Get alive neighbour vertices of `v` (sorted, deduped).
     *
     *  Production code path uses `try_collapse`'s `outKeepNbrs` out-param to
     *  avoid this method's per-call allocation; this overload is retained
     *  only as the oracle for the `HLSCMInternal.TryCollapse_*` tests.
     */
    [[nodiscard]] auto vertex_neighbors(std::size_t v) const -> std::vector<std::size_t>
    {
        std::vector<std::size_t> nbrs;
        for (auto fi : vert_faces_[v]) {
            if (!face_alive_[fi])
                continue;
            for (auto vi : faces_[fi]) {
                if (vi != v && alive_[vi])
                    nbrs.push_back(vi);
            }
        }
        std::sort(nbrs.begin(), nbrs.end());
        nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
        return nbrs;
    }

    /** Take a snapshot of surviving vertices and faces for a hierarchy level */
    [[nodiscard]] auto snapshot() const -> HierarchyLevel<T>
    {
        HierarchyLevel<T> level;
        level.original_to_local.assign(alive_.size(), std::nullopt);

        // Build mapping from original indices to level-local indices
        std::size_t localIdx = 0;
        for (std::size_t i = 0; i < alive_.size(); ++i) {
            if (alive_[i]) {
                level.original_to_local[i] = localIdx;
                level.local_to_original.push_back(i);
                level.positions.push_back(positions_[i]);
                localIdx++;
            }
        }

        // Remap faces. Every alive face references only alive vertices, so
        // the unwrap is always safe here.
        for (std::size_t fi = 0; fi < faces_.size(); ++fi) {
            if (!face_alive_[fi]) {
                continue;
            }
            std::array<std::size_t, 3> localTri;
            for (int j = 0; j < 3; ++j) {
                localTri[j] = *level.original_to_local[faces_[fi][j]];
            }
            level.faces.push_back(localTri);
        }

        return level;
    }

    /** Rebuild the edge list from alive faces and return it */
    auto rebuild_and_get_edges() -> const std::vector<std::pair<std::size_t, std::size_t>>&
    {
        build_edges_();
        return edges_;
    }

private:
    /** Recompute per-vertex QEM quadrics from the live faces (Garland-Heckbert). */
    void compute_quadrics_()
    {
        for (auto& q : quadrics_) {
            q = Quadric<T>();
        }

        for (std::size_t fi = 0; fi < faces_.size(); ++fi) {
            if (!face_alive_[fi]) {
                continue;
            }
            auto& tri = faces_[fi];
            auto& p0 = positions_[tri[0]];
            auto& p1 = positions_[tri[1]];
            auto& p2 = positions_[tri[2]];

            // Face plane: normal = (p1-p0) x (p2-p0), normalized
            auto e1 = p1 - p0;
            auto e2 = p2 - p0;
            auto n = cross(e1, e2);
            auto len = norm(n);
            if (len < std::numeric_limits<T>::epsilon()) {
                continue;
            }
            n /= len;

            T a = n[0], b = n[1], c = n[2];
            T d = -(a * p0[0] + b * p0[1] + c * p0[2]);

            Quadric<T> faceQ(a, b, c, d);
            quadrics_[tri[0]] += faceQ;
            quadrics_[tri[1]] += faceQ;
            quadrics_[tri[2]] += faceQ;
        }
    }

    /** Rebuild the unique-edge list from the live faces. */
    void build_edges_()
    {
        edges_.clear();
        std::unordered_set<std::uint64_t> seen;
        auto edgeKey = [this](std::size_t a, std::size_t b) -> std::uint64_t {
            auto n = static_cast<std::uint64_t>(positions_.size());
            return static_cast<std::uint64_t>(std::min(a, b)) * n +
                   static_cast<std::uint64_t>(std::max(a, b));
        };

        for (std::size_t fi = 0; fi < faces_.size(); ++fi) {
            if (!face_alive_[fi]) {
                continue;
            }
            auto& tri = faces_[fi];
            for (int j = 0; j < 3; ++j) {
                auto a = tri[j];
                auto b = tri[(j + 1) % 3];
                auto key = edgeKey(a, b);
                if (seen.insert(key).second) {
                    edges_.emplace_back(std::min(a, b), std::max(a, b));
                }
            }
        }
    }

    /** Compute barycentric coordinates of point p in triangle (a, b, c) */
    static auto computeBarycentric_(const Vec<T, 3>& a, const Vec<T, 3>& b, const Vec<T, 3>& c,
                                    const Vec<T, 3>& p) -> std::array<T, 3>
    {
        auto v0 = b - a;
        auto v1 = c - a;
        auto v2 = p - a;

        T d00 = dot(v0, v0);
        T d01 = dot(v0, v1);
        T d11 = dot(v1, v1);
        T d20 = dot(v2, v0);
        T d21 = dot(v2, v1);

        T denom = d00 * d11 - d01 * d01;
        if (std::abs(denom) < std::numeric_limits<T>::epsilon()) {
            return {T(1), T(0), T(0)};
        }

        T v = (d11 * d20 - d01 * d21) / denom;
        T w = (d00 * d21 - d01 * d20) / denom;

        // Clamp to avoid extreme extrapolation from far-away collapses
        v = std::max(T(-0.5), std::min(T(1.5), v));
        w = std::max(T(-0.5), std::min(T(1.5), w));
        T u = T(1) - v - w;

        return {u, v, w};
    }

    /** Per-vertex 3D positions, indexed by original vertex index. */
    std::vector<Vec<T, 3>> positions_;
    /** Per-vertex alive flag; collapses set entries to false. */
    std::vector<bool> alive_;
    /** Per-vertex boundary flag, copied from the source mesh. */
    std::vector<bool> is_boundary_;
    /** Per-vertex pinned flag; pinned vertices cannot be removed by collapse. */
    std::vector<bool> is_pinned_;
    /** Per-vertex QEM quadrics accumulated from incident face planes. */
    std::vector<Quadric<T>> quadrics_;
    /** Face connectivity: each face is three vertex indices into `positions_`. */
    std::vector<std::array<std::size_t, 3>> faces_;
    /** Per-face alive flag; collapses retire faces by setting entries to false. */
    std::vector<bool> face_alive_;
    /** Per-vertex incident face list. */
    std::vector<std::vector<std::size_t>> vert_faces_;
    /** Count of alive vertices (kept current across collapses). */
    std::size_t num_alive_verts_{0};
    /** Count of alive faces (kept current across collapses). */
    std::size_t num_alive_faces_{0};
    /** Unique-edge list, rebuilt on demand by `build_edges_()`. */
    std::vector<std::pair<std::size_t, std::size_t>> edges_;

    /** Scratch buffer: neighbor indices shared between two endpoints. */
    std::vector<std::size_t> scratch_shared_;
    /** Scratch buffer: faces to retire on a collapse. */
    std::vector<std::size_t> scratch_remove_;
    /** Scratch buffer: post-collapse face rewrites. */
    std::vector<std::array<std::size_t, 3>> scratch_post_faces_;
    /** Scratch buffer: neighbor list of one endpoint. */
    std::vector<std::size_t> scratch_nbrs_a_;
    /** Scratch buffer: neighbor list of the other endpoint. */
    std::vector<std::size_t> scratch_nbrs_b_;
    /** Scratch buffer: deduplicated union of `scratch_nbrs_a_` and `scratch_nbrs_b_`. */
    std::vector<std::size_t> scratch_shared_nbrs_;
};

/**
 * @brief Build a mesh hierarchy by greedy QEM decimation
 *
 * Returns a vector of HierarchyLevel from finest to coarsest, plus
 * the collapse records needed for prolongation (ordered from finest to coarsest).
 */
template <typename T, class MeshPtr>
auto BuildHierarchy(const MeshPtr& mesh, const std::vector<std::size_t>& pinIndices,
                    std::size_t levelRatio, std::size_t minCoarseVerts)
    -> std::pair<std::vector<HierarchyLevel<T>>, std::vector<std::vector<CollapseRecord<T>>>>
{
    DecimationMesh<T> dmesh;
    dmesh.build(mesh, pinIndices);

    // Finest level snapshot
    std::vector<HierarchyLevel<T>> levels;
    levels.push_back(dmesh.snapshot());

    std::vector<std::vector<CollapseRecord<T>>> collapsesByLevel;

    auto targetVerts = dmesh.num_alive_verts();
    if (targetVerts <= minCoarseVerts) {
        // Mesh is already small enough — single level
        return {levels, collapsesByLevel};
    }

    // Reused across collapses to avoid per-collapse heap allocation; populated
    // by try_collapse with vKeep's post-collapse neighbors.
    std::vector<std::size_t> nbrs;

    while (targetVerts > minCoarseVerts) {
        auto nextTarget = std::max(targetVerts / levelRatio, minCoarseVerts);
        std::vector<CollapseRecord<T>> levelCollapses;

        // Build priority queue of edge collapses
        using CostEdge = std::pair<T, std::pair<std::size_t, std::size_t>>;
        std::priority_queue<CostEdge, std::vector<CostEdge>, std::greater<CostEdge>> pq;

        const auto& edges = dmesh.rebuild_and_get_edges();
        for (auto& [a, b] : edges) {
            // Try collapsing the collapsible vertex towards the other
            if (dmesh.is_collapsible(a)) {
                pq.push({dmesh.collapse_cost(a, b), {a, b}});
            }
            if (dmesh.is_collapsible(b)) {
                pq.push({dmesh.collapse_cost(b, a), {b, a}});
            }
        }

        while (dmesh.num_alive_verts() > nextTarget && !pq.empty()) {
            auto [cost, edge] = pq.top();
            pq.pop();

            auto [vRemove, vKeep] = edge;
            // Skip stale entries whose endpoints have already been collapsed
            if (!dmesh.is_alive(vRemove) || !dmesh.is_alive(vKeep)) {
                continue;
            }
            auto record = dmesh.try_collapse(vRemove, vKeep, &nbrs);
            if (!record) {
                continue;
            }

            levelCollapses.push_back(*record);

            // Add new edges involving vKeep to the priority queue using the
            // post-collapse neighbor list populated by try_collapse.
            for (auto nb : nbrs) {
                if (dmesh.is_collapsible(nb)) {
                    pq.push({dmesh.collapse_cost(nb, vKeep), {nb, vKeep}});
                }
                if (dmesh.is_collapsible(vKeep)) {
                    pq.push({dmesh.collapse_cost(vKeep, nb), {vKeep, nb}});
                }
            }
        }

        if (levelCollapses.empty()) {
            break;  // No more valid collapses possible
        }

        collapsesByLevel.push_back(std::move(levelCollapses));
        levels.push_back(dmesh.snapshot());
        targetVerts = dmesh.num_alive_verts();
    }

    return {levels, collapsesByLevel};
}

/**
 * @brief Build a HalfEdgeMesh from a hierarchy level
 */
template <typename T>
auto BuildLevelMesh(const HierarchyLevel<T>& level) -> typename HalfEdgeMesh<T>::Pointer
{
    auto mesh = HalfEdgeMesh<T>::New();

    for (const auto& pos : level.positions) {
        mesh->insert_vertex(pos[0], pos[1], pos[2]);
    }

    // level.faces is already `vector<array<size_t,3>>`; insert_faces is generic
    // over containers-of-iterables so we can pass it directly instead of
    // copying every face into a `vector<vector<size_t>>` (one heap allocation
    // per face).
    mesh->insert_faces(level.faces);

    return mesh;
}

/** UV vector indexed by original (finest-level) vertex idx. Slots for
 *  vertices not yet solved/prolongated hold `std::nullopt`. */
template <typename T>
using UVVector = std::vector<std::optional<Vec<T, 2>>>;

/**
 * @brief Prolongate UV coordinates from a coarser level to a finer level
 *
 * Surviving vertices keep their existing UVs; vertices that were removed by
 * a collapse get UVs via barycentric interpolation in their containing
 * post-collapse triangle.
 *
 * @param uvs UV coordinates indexed by original vertex index. Unset slots are
 *            `std::nullopt`. Mutated in place: every vertex removed in
 *            `collapses` is filled in. Returned by move.
 * @param collapses Collapse records for this level transition.
 */
template <typename T>
auto ProlongateUVs(UVVector<T> uvs, const std::vector<CollapseRecord<T>>& collapses) -> UVVector<T>
{
    // Undo collapses in reverse order — the last collapse applied is the first
    // we need to undo to recover the next-finer level's UVs. All three
    // containing-tri vertices are guaranteed to have UVs by the time we
    // dereference them: they survived the collapse we're undoing.
    for (auto it = collapses.rbegin(); it != collapses.rend(); ++it) {
        auto& rec = *it;
        auto& tri = rec.containing_tri;

        uvs[rec.v_removed] =
            *uvs[tri[0]] * rec.bary[0] + *uvs[tri[1]] * rec.bary[1] + *uvs[tri[2]] * rec.bary[2];
    }

    return uvs;
}

/**
 * @brief Solve the LSCM system at one hierarchy level
 *
 * Builds the Lévy et al. Eq. 10 LSCM system on the given level mesh with the
 * given pin vertices. If an initial guess is provided, uses solveWithGuess.
 *
 * @param origVertCount Original (finest-level) vertex count; sizes the output
 *                      UV vector so it can be indexed by original vertex idx.
 * @return UV coordinates indexed by original vertex index. Slots for vertices
 *         not present at this level remain `std::nullopt`.
 */
template <typename T, class SolverType>
auto SolveLSCMLevel(const typename HalfEdgeMesh<T>::Pointer& levelMesh,
                    const detail::hlscm::HierarchyLevel<T>& level,
                    const detail::lscm::PinMap<T>& origPins, std::size_t origVertCount,
                    const UVVector<T>* initialGuess) -> UVVector<T>
{
    using SparseMatrix = Eigen::SparseMatrix<T>;
    using DenseMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using Mesh = HalfEdgeMesh<T>;

    auto numVerts = levelMesh->num_vertices();
    auto numFixed = origPins.size();
    auto numFree = numVerts - numFixed;

    // Map each original pin idx to its level-local idx. Pins are guaranteed
    // to survive every decimation level (DecimationMesh::try_collapse rejects
    // collapses of pinned vertices), so this unwrap is expected to succeed.
    // Throw a meaningful exception rather than std::bad_optional_access if a
    // future refactor ever breaks the invariant.
    detail::lscm::PinMap<T> localPins;
    localPins.reserve(numFixed);
    std::unordered_set<std::size_t> localPinIdx;
    localPinIdx.reserve(numFixed);
    for (const auto& [origIdx, uv] : origPins) {
        const auto& localOpt = level.original_to_local[origIdx];
        if (!localOpt.has_value()) {
            throw SolverException(
                "HLSCM: pinned vertex was removed during decimation (invariant violated)");
        }
        auto localIdx = *localOpt;
        localPins.emplace_back(localIdx, uv);
        localPinIdx.insert(localIdx);
    }

    // Pin placement + LSCM system assembly (shared with AngleBasedLSCM)
    auto parts = detail::lscm::BuildSystem<T, Mesh>(levelMesh, localPins);
    auto& A = parts.A;
    auto& b = parts.b;
    auto& freeIdxTable = parts.freeIdxTable;

    // Build initial guess vector from prolongated UVs.
    // `initialGuess` is indexed by original vertex idx; nullopt entries are
    // vertices not yet solved at any coarser level.
    bool warmed = initialGuess != nullptr;
    auto buildInitialGuess = [&]() -> DenseMatrix {
        DenseMatrix x0 = DenseMatrix::Zero(2 * numFree, 1);
        for (const auto& v : levelMesh->vertices()) {
            if (localPinIdx.count(v->idx)) {
                continue;
            }
            auto freeIdx = freeIdxTable.at(v->idx);
            auto origIdx = level.local_to_original[v->idx];
            if (const auto& guess = (*initialGuess)[origIdx]) {
                x0(2 * freeIdx, 0) = (*guess)[0];
                x0(2 * freeIdx + 1, 0) = (*guess)[1];
            }
        }
        return x0;
    };

    // Convergence tolerance for iterative solvers.  Eigen defaults to
    // machine epsilon (~2e-16 for double) which is far tighter than
    // needed for UV parameterization and prevents the warm-start from
    // reducing iteration count.  1e-8 gives ~8 digits of relative
    // residual precision — more than sufficient for texturing.
    constexpr T kTolerance = T(1e-8);

    // Solve
    DenseMatrix x;
    if constexpr (detail::is_instance_of_v<SolverType, Eigen::LeastSquaresConjugateGradient>) {
        // LSCG operates on the rectangular system A directly (avoids squaring the condition number)
        SolverType solver(A);
        solver.setTolerance(kTolerance);
        DenseMatrix bDense = DenseMatrix(b);
        if (warmed) {
            x = solver.solveWithGuess(bDense, buildInitialGuess());
        } else {
            x = solver.solve(bDense);
        }
        if (solver.info() == Eigen::ComputationInfo::NumericalIssue ||
            solver.info() == Eigen::ComputationInfo::InvalidInput ||
            solver.info() == Eigen::ComputationInfo::NoConvergence) {
            throw SolverException("HLSCM: LSCG solve failed at hierarchy level");
        }
    } else if constexpr (std::is_base_of_v<Eigen::IterativeSolverBase<SolverType>, SolverType>) {
        // CG and other iterative solvers on the square SPD system AtA.
        SparseMatrix AtA = A.transpose() * A;
        AtA.makeCompressed();
        DenseMatrix Atb = DenseMatrix(A.transpose() * b);
        SolverType solver(AtA);
        solver.setTolerance(kTolerance);
        if (warmed) {
            x = solver.solveWithGuess(Atb, buildInitialGuess());
        } else {
            x = solver.solve(Atb);
        }
        if (solver.info() == Eigen::ComputationInfo::NumericalIssue ||
            solver.info() == Eigen::ComputationInfo::InvalidInput ||
            solver.info() == Eigen::ComputationInfo::NoConvergence) {
            throw SolverException("HLSCM: iterative solve failed at hierarchy level");
        }
    } else {
        // Direct solver: decompose AtA and solve. No warm-start benefit.
        SparseMatrix AtA = A.transpose() * A;
        AtA.makeCompressed();
        DenseMatrix Atb = DenseMatrix(A.transpose() * b);
        SolverType solver;
        solver.compute(AtA);
        if (solver.info() != Eigen::ComputationInfo::Success) {
            throw SolverException("HLSCM: solver decomposition failed");
        }
        x = solver.solve(Atb);
    }

    // Build output UV vector indexed by original vertex idx. Vertices not
    // present at this level remain nullopt; they will be filled in by
    // prolongation when undoing collapses at finer levels.
    UVVector<T> uvs(origVertCount);
    for (const auto& [origIdx, uv] : origPins) {
        uvs[origIdx] = uv;
    }
    for (const auto& v : levelMesh->vertices()) {
        if (localPinIdx.count(v->idx)) {
            continue;
        }
        auto freeIdx = 2 * freeIdxTable.at(v->idx);
        auto origIdx = level.local_to_original[v->idx];
        uvs[origIdx] = Vec<T, 2>(x(freeIdx, 0), x(freeIdx + 1, 0));
    }
    return uvs;
}

}  // namespace hlscm
}  // namespace detail

/**
 * @brief Compute parameterized mesh using Hierarchical LSCM
 *
 * Implements the HLSCM algorithm from Ray & Lévy, "Hierarchical Least Squares
 * Conformal Map" (2003) \cite ray2003hlscm. Uses cascadic multigrid to
 * accelerate LSCM
 * convergence: the mesh is decimated into a hierarchy, LSCM is solved on the
 * coarsest level, and the solution is prolongated and refined at each finer
 * level using conjugate gradient with the prolongated UVs as initial guess.
 *
 * For small meshes the hierarchy has a single level and HLSCM degrades
 * gracefully to a standard LSCM solve.
 *
 * @tparam T Floating-point type
 * @tparam MeshType HalfEdgeMesh type which implements the default mesh traits
 * @tparam Solver An Eigen iterative or direct solver. Iterative solvers
 *         (ConjugateGradient, LeastSquaresConjugateGradient) support warm-
 *         starting from the coarser-level solution; direct solvers ignore the
 *         initial guess. Defaults to
 *         `ConjugateGradient<SparseMatrix<T>, Lower|Upper>`, which solves the
 *         normal equations (AᵀA x = Aᵀb) of the overdetermined LSCM system.
 *         The `Lower|Upper` flag enables OpenMP-parallelized SpMV on the
 *         symmetric AᵀA matrix, giving the best multi-thread performance.
 *         LeastSquaresConjugateGradient is a valid alternative; it operates
 *         on the same normal equations internally but without the OpenMP
 *         benefit.
 *
 * @note **Solver default differs from AngleBasedLSCM.** AngleBasedLSCM
 *       defaults to SparseLU (a direct solver); HierarchicalLSCM defaults
 *       to ConjugateGradient so it can warm-start from the coarser-level
 *       solution. Using a direct solver via the `Solver` template parameter
 *       is valid but disables warm-starting — the initial guess is ignored.
 *
 * @note **Single-level fallback.** When the mesh is too small to decimate
 *       (all vertices are boundary, or minCoarseVertices is already reached),
 *       HLSCM falls back to a standard LSCM solve using *this class's*
 *       `Solver` template parameter, not AngleBasedLSCM's default (SparseLU).
 *       The result is numerically equivalent but may differ in convergence
 *       behavior from a plain AngleBasedLSCM call.
 */
template <typename T, class MeshType = HalfEdgeMesh<T>,
          class Solver =
              Eigen::ConjugateGradient<Eigen::SparseMatrix<T>, Eigen::Lower | Eigen::Upper>,
          std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
class HierarchicalLSCM
{
public:
    /** @brief Mesh type alias */
    using Mesh = MeshType;

    /**
     * @brief Per-pin entry: (mesh vertex index, target UV).
     *
     * A PinMap of size N ≥ 2 specifies an explicit LSCM pin set. Each pin
     * survives every decimation level — `DecimationMesh::try_collapse` rejects
     * any collapse that would remove a pinned vertex.
     */
    using PinMap = detail::lscm::PinMap<T>;

    /**
     * @brief Set the explicit pin set used by `compute()`
     */
    void set_pins(PinMap pins)
    {
        pins_ = std::move(pins);
        legacy_pin_indices_.reset();
    }

    /**
     * @brief Deprecated: set a pin pair by index using the LSCM axis-snap
     * convention at compute time.
     *
     * @deprecated Prefer `set_pins(PinMap)`. This overload will be removed in
     * version 3.0.
     */
    [[deprecated("Use set_pins(PinMap); will be removed in 3.0")]] void setPinnedVertices(
        std::size_t pin0Idx, std::size_t pin1Idx)
    {
        legacy_pin_indices_ = {pin0Idx, pin1Idx};
        pins_.reset();
    }

    /** @brief Set the vertex ratio between consecutive hierarchy levels (default: 10) */
    void set_level_ratio(std::size_t ratio)
    {
        if (ratio < 2) {
            throw std::invalid_argument("HierarchicalLSCM: level_ratio must be >= 2");
        }
        level_ratio_ = ratio;
    }

    /** @brief Set the minimum vertex count at the coarsest level (default: 100) */
    void set_min_coarse_vertices(std::size_t count)
    {
        if (count < 3) {
            throw std::invalid_argument("HierarchicalLSCM: min_coarse_vertices must be >= 3");
        }
        min_coarse_vertices_ = count;
    }

    /** @deprecated Use `set_level_ratio`; will be removed in 3.0. */
    [[deprecated("Use set_level_ratio; will be removed in 3.0")]] void setLevelRatio(
        std::size_t ratio)
    {
        set_level_ratio(ratio);
    }

    /** @deprecated Use `set_min_coarse_vertices`; will be removed in 3.0. */
    [[deprecated("Use set_min_coarse_vertices; will be removed in 3.0")]] void setMinCoarseVertices(
        std::size_t count)
    {
        set_min_coarse_vertices(count);
    }

    /**
     * @brief Compute parameterization using instance configuration
     *
     * If `set_pins()` was called, uses the supplied PinMap. Otherwise, if the
     * deprecated `setPinnedVertices()` was called, builds a PinMap from those
     * indices via the LSCM axis-snap convention. Otherwise auto-selects two
     * boundary vertices using the same logic as `Compute(mesh)`.
     *
     * @throws MeshException if pin selection fails (no boundary vertices)
     * @throws SolverException if any hierarchy level fails to solve
     */
    void compute(typename Mesh::Pointer& mesh) const
    {
        PinMap pins;
        if (pins_) {
            pins = *pins_;
            detail::lscm::ValidatePins<T, Mesh>(mesh, pins);
        } else if (legacy_pin_indices_) {
            pins = detail::lscm::AutoPlacePair<T, Mesh>(mesh, legacy_pin_indices_->first,
                                                        legacy_pin_indices_->second);
        } else {
            pins = detail::lscm::AutoSelectPins<T, Mesh>(mesh);
        }
        ComputeImpl(mesh, pins, level_ratio_, min_coarse_vertices_);
    }

    /**
     * @brief Compute with automatic pin selection
     *
     * Selects pins identically to AngleBasedLSCM::Compute() — first boundary
     * vertex and its boundary-edge neighbor, placed via LSCM axis-snap.
     *
     * @throws MeshException if the mesh has no boundary vertices (pin selection
     *         fails) or the mesh is otherwise invalid
     * @throws SolverException if any hierarchy level fails to solve
     */
    static void Compute(typename Mesh::Pointer& mesh)
    {
        ComputeImpl(mesh, detail::lscm::AutoSelectPins<T, Mesh>(mesh));
    }

    /**
     * @brief Compute with caller-specified pin UVs
     *
     * @throws std::invalid_argument If `pins` has fewer than two entries, a
     * duplicate index, or an out-of-range index.
     * @throws SolverException If any hierarchy level fails to solve.
     */
    static void Compute(typename Mesh::Pointer& mesh, const PinMap& pins)
    {
        detail::lscm::ValidatePins<T, Mesh>(mesh, pins);
        ComputeImpl(mesh, pins);
    }

    /**
     * @brief Deprecated: compute with an explicit pin pair by index.
     *
     * Builds a 2-entry PinMap using the LSCM axis-snap convention and
     * dispatches to the PinMap path.
     *
     * @deprecated Prefer `Compute(mesh, PinMap)`. This overload will be
     * removed in version 3.0.
     */
    [[deprecated("Use Compute(mesh, PinMap); will be removed in 3.0")]] static void Compute(
        typename Mesh::Pointer& mesh, std::size_t pin0Idx, std::size_t pin1Idx)
    {
        ComputeImpl(mesh, detail::lscm::AutoPlacePair<T, Mesh>(mesh, pin0Idx, pin1Idx));
    }

private:
    /**
     * @brief Copy edge angles from the original mesh to a level mesh
     *
     * At the finest hierarchy level (k=0), the level mesh has the same
     * face/edge structure as the original mesh. If ABF was run beforehand,
     * we must use the ABF-optimized angles rather than recomputing from
     * geometry. Coarser levels always use geometry angles.
     */
    static void CopyAnglesFromOriginal(const typename Mesh::Pointer& original,
                                       const typename HalfEdgeMesh<T>::Pointer& levelMesh)
    {
        for (const auto& f : levelMesh->faces()) {
            auto origFace = original->face(f->idx);
            auto le = f->head;
            auto oe = origFace->head;
            for (int j = 0; j < 3; ++j) {
                le->alpha = oe->alpha;
                le = le->next;
                oe = oe->next;
            }
        }
    }

    /**
     * @brief Core hierarchical LSCM solve given a resolved PinMap
     *
     * Builds the mesh hierarchy with each pinned vertex flagged as
     * non-collapsible, solves LSCM at the coarsest level, then prolongates and
     * refines at each finer level. Pin UVs come from the PinMap verbatim;
     * non-pin vertex UVs are written by the solve.
     */
    static void ComputeImpl(typename Mesh::Pointer& mesh, const PinMap& pins,
                            std::size_t levelRatio = 10, std::size_t minCoarseVerts = 100)
    {
        // Extract just the indices for the hierarchy's non-collapsible set.
        std::vector<std::size_t> pinIndices;
        pinIndices.reserve(pins.size());
        for (const auto& [vIdx, uv] : pins) {
            pinIndices.push_back(vIdx);
        }

        // Build mesh hierarchy
        auto [levels, collapsesByLevel] =
            detail::hlscm::BuildHierarchy<T>(mesh, pinIndices, levelRatio, minCoarseVerts);

        if (levels.size() <= 1) {
            // Mesh too small for hierarchy — single-level LSCM solve.
            // Delegate to AngleBasedLSCM with the same pins.
            AngleBasedLSCM<T, MeshType, Solver>::Compute(mesh, pins);
            return;
        }

        // UV vectors are sized to the original (finest-level) vertex count
        // so they can be indexed by `original vertex idx` at every level.
        const auto origVertCount = mesh->num_vertices();

        // Solve coarsest level (last in the array)
        auto coarsestIdx = levels.size() - 1;
        auto coarseMesh = detail::hlscm::BuildLevelMesh<T>(levels[coarsestIdx]);
        ComputeMeshAngles(coarseMesh);
        auto uvs = detail::hlscm::SolveLSCMLevel<T, Solver>(coarseMesh, levels[coarsestIdx], pins,
                                                            origVertCount, nullptr);

        // Prolongate and refine at each finer level
        for (std::size_t k = coarsestIdx; k-- > 0;) {
            // Prolongate UVs from level k+1 to level k
            uvs = detail::hlscm::ProlongateUVs<T>(std::move(uvs), collapsesByLevel[k]);

            // Build level mesh
            auto levelMesh = detail::hlscm::BuildLevelMesh<T>(levels[k]);

            if (k == 0) {
                // Finest level: use original mesh angles (may be ABF-optimized)
                CopyAnglesFromOriginal(mesh, levelMesh);
            } else {
                // Coarser levels: compute angles from 3D geometry
                ComputeMeshAngles(levelMesh);
            }

            // Solve with initial guess
            uvs = detail::hlscm::SolveLSCMLevel<T, Solver>(levelMesh, levels[k], pins,
                                                           origVertCount, &uvs);
        }

        // Transfer final UVs back to input mesh.
        // At the finest level every vertex has a UV; the unwrap is safe.
        for (const auto& v : mesh->vertices()) {
            const auto& uv = *uvs[v->idx];
            v->pos[0] = uv[0];
            v->pos[1] = uv[1];
            v->pos[2] = T(0);
        }
    }

    /** Optional explicit pin set configured via `set_pins()`. */
    std::optional<PinMap> pins_;
    /** Deprecated: legacy two-pin index pair set via `setPinnedVertices`. */
    std::optional<std::pair<std::size_t, std::size_t>> legacy_pin_indices_;
    /** Ratio of vertices between consecutive hierarchy levels */
    std::size_t level_ratio_{10};
    /** Minimum vertex count for the coarsest hierarchy level */
    std::size_t min_coarse_vertices_{100};
};

}  // namespace OpenABF
