#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <optional>
#include <queue>
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
namespace hlscm
{

/** @brief Symmetric 4×4 quadric matrix for QEM error metric (Garland-Heckbert) */
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

    auto operator+=(const Quadric& o) -> Quadric&
    {
        for (std::size_t i = 0; i < 10; ++i) {
            q[i] += o.q[i];
        }
        return *this;
    }

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

/** @brief Record of a single half-edge collapse for prolongation */
template <typename T>
struct CollapseRecord {
    /** Index of the removed vertex (in the original/fine mesh) */
    std::size_t vRemoved;
    /** Index of the kept vertex (in the original/fine mesh) */
    std::size_t vKept;
    /** Post-collapse triangle containing vRemoved (original vertex indices) */
    std::array<std::size_t, 3> containingTri;
    /** Barycentric coordinates of vRemoved in containingTri */
    std::array<T, 3> bary;
};

/** @brief A level in the mesh hierarchy */
template <typename T>
struct HierarchyLevel {
    /** Vertex positions (indexed by level-local index) */
    std::vector<Vec<T, 3>> positions;
    /** Face connectivity (each face is 3 level-local indices) */
    std::vector<std::array<std::size_t, 3>> faces;
    /** Map from level-local vertex index to original (finest) vertex index */
    std::vector<std::size_t> localToOriginal;
    /** Map from original vertex index to level-local index */
    std::unordered_map<std::size_t, std::size_t> originalToLocal;
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
    void build(const MeshPtr& mesh, std::size_t pin0, std::size_t pin1)
    {
        auto nv = mesh->num_vertices();
        auto nf = mesh->num_faces();

        positions_.resize(nv);
        alive_.assign(nv, true);
        isBoundary_.assign(nv, false);
        isPinned_.assign(nv, false);
        quadrics_.resize(nv);

        isPinned_[pin0] = true;
        isPinned_[pin1] = true;

        for (const auto& v : mesh->vertices()) {
            positions_[v->idx] = v->pos;
            isBoundary_[v->idx] = v->is_boundary();
        }

        faces_.reserve(nf);
        faceAlive_.reserve(nf);
        vertFaces_.resize(nv);

        for (const auto& f : mesh->faces()) {
            auto e0 = f->head;
            auto e1 = e0->next;
            auto e2 = e1->next;
            std::array<std::size_t, 3> tri{e0->vertex->idx, e1->vertex->idx, e2->vertex->idx};
            auto fi = faces_.size();
            faces_.push_back(tri);
            faceAlive_.push_back(true);
            vertFaces_[tri[0]].push_back(fi);
            vertFaces_[tri[1]].push_back(fi);
            vertFaces_[tri[2]].push_back(fi);
        }

        numAliveVerts_ = nv;
        numAliveFaces_ = nf;

        computeQuadrics_();
        buildEdges_();
    }

    /** Get number of alive vertices */
    [[nodiscard]] auto numAliveVerts() const -> std::size_t { return numAliveVerts_; }

    /** Get number of alive faces */
    [[nodiscard]] auto numAliveFaces() const -> std::size_t { return numAliveFaces_; }

    /**
     * @brief Try to collapse edge (vRemove → vKeep), returning a collapse record
     *
     * Returns nullopt if the collapse is invalid.
     */
    auto tryCollapse(std::size_t vRemove, std::size_t vKeep) -> std::optional<CollapseRecord<T>>
    {
        if (!alive_[vRemove] || !alive_[vKeep]) {
            return std::nullopt;
        }
        if (isBoundary_[vRemove] || isPinned_[vRemove]) {
            return std::nullopt;
        }

        // Find shared faces (will be removed) and vRemove-only faces (will be updated)
        std::vector<std::size_t> sharedFaces;
        std::vector<std::size_t> removeFaces;
        for (auto fi : vertFaces_[vRemove]) {
            if (!faceAlive_[fi]) {
                continue;
            }
            bool hasKeep = false;
            for (auto vi : faces_[fi]) {
                if (vi == vKeep) {
                    hasKeep = true;
                    break;
                }
            }
            if (hasKeep) {
                sharedFaces.push_back(fi);
            } else {
                removeFaces.push_back(fi);
            }
        }

        // Exactly 2 shared faces for interior edge in manifold mesh
        // (could be 1 for boundary, but we don't collapse boundary vertices)
        if (sharedFaces.size() != 2) {
            return std::nullopt;
        }

        // Check for topology validity: the link condition
        // Vertices shared between vRemove's and vKeep's neighborhoods
        // (excluding vRemove and vKeep themselves) must be exactly the
        // two vertices opposite the shared faces.
        auto neighborsOf = [&](std::size_t v) {
            std::unordered_set<std::size_t> nbrs;
            for (auto fi : vertFaces_[v]) {
                if (!faceAlive_[fi]) {
                    continue;
                }
                for (auto vi : faces_[fi]) {
                    if (vi != v) {
                        nbrs.insert(vi);
                    }
                }
            }
            return nbrs;
        };

        auto nbrsRemove = neighborsOf(vRemove);
        auto nbrsKeep = neighborsOf(vKeep);

        std::unordered_set<std::size_t> sharedNbrs;
        for (auto v : nbrsRemove) {
            if (v != vKeep && nbrsKeep.count(v)) {
                sharedNbrs.insert(v);
            }
        }

        // The shared neighbors should be exactly the "opposite" vertices of
        // the shared faces
        std::unordered_set<std::size_t> expectedShared;
        for (auto fi : sharedFaces) {
            for (auto vi : faces_[fi]) {
                if (vi != vRemove && vi != vKeep) {
                    expectedShared.insert(vi);
                }
            }
        }
        if (sharedNbrs != expectedShared) {
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
        std::vector<std::array<std::size_t, 3>> postFaces;
        // Existing vKeep faces (excluding shared faces which will be removed)
        std::unordered_set<std::size_t> sharedSet(sharedFaces.begin(), sharedFaces.end());
        for (auto fi : vertFaces_[vKeep]) {
            if (!faceAlive_[fi] || sharedSet.count(fi)) {
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
        // Store bary coords of vRemoved's position in the post-collapse triangle.
        CollapseRecord<T> record;
        record.vRemoved = vRemove;
        record.vKept = vKeep;

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
            record.containingTri = postTri;
            record.bary = computeBarycentric_(positions_[postTri[0]], positions_[postTri[1]],
                                              positions_[postTri[2]], positions_[vRemove]);
        } else {
            // Edge case: all faces are shared — vertex collapses directly onto vKeep
            record.containingTri = {vKeep, vKeep, vKeep};
            record.bary = {T(1), T(0), T(0)};
        }

        // Execute collapse
        // Kill shared faces
        for (auto fi : sharedFaces) {
            faceAlive_[fi] = false;
            numAliveFaces_--;
        }

        // Update removeFaces: replace vRemove with vKeep
        for (auto fi : removeFaces) {
            for (auto& vi : faces_[fi]) {
                if (vi == vRemove) {
                    vi = vKeep;
                }
            }
            vertFaces_[vKeep].push_back(fi);
        }

        // Mark vRemove as dead
        alive_[vRemove] = false;
        numAliveVerts_--;

        // Merge quadrics
        quadrics_[vKeep] += quadrics_[vRemove];

        // Compact dead face indices from vKeep's adjacency list
        auto& vkFaces = vertFaces_[vKeep];
        vkFaces.erase(std::remove_if(vkFaces.begin(), vkFaces.end(),
                                     [this](std::size_t fi) { return !faceAlive_[fi]; }),
                      vkFaces.end());

        return record;
    }

    /** Compute collapse cost for edge (v0 → v1): Q_merged evaluated at v1 */
    [[nodiscard]] auto collapseCost(std::size_t v0, std::size_t v1) const -> T
    {
        auto Q = quadrics_[v0] + quadrics_[v1];
        auto& p = positions_[v1];
        return Q.evaluate(p[0], p[1], p[2]);
    }

    /** Check if a vertex is alive */
    [[nodiscard]] auto isAlive(std::size_t v) const -> bool { return alive_[v]; }

    /** Check if a vertex is collapsible (not boundary, not pinned, alive) */
    [[nodiscard]] auto isCollapsible(std::size_t v) const -> bool
    {
        return alive_[v] && !isBoundary_[v] && !isPinned_[v];
    }

    /** Get edges incident to vertex v (pairs of (v, neighbor)) */
    [[nodiscard]] auto vertexNeighbors(std::size_t v) const -> std::vector<std::size_t>
    {
        std::unordered_set<std::size_t> nbrs;
        for (auto fi : vertFaces_[v]) {
            if (!faceAlive_[fi]) {
                continue;
            }
            for (auto vi : faces_[fi]) {
                if (vi != v && alive_[vi]) {
                    nbrs.insert(vi);
                }
            }
        }
        return {nbrs.begin(), nbrs.end()};
    }

    /** Take a snapshot of surviving vertices and faces for a hierarchy level */
    [[nodiscard]] auto snapshot() const -> HierarchyLevel<T>
    {
        HierarchyLevel<T> level;

        // Build mapping from original indices to level-local indices
        std::size_t localIdx = 0;
        for (std::size_t i = 0; i < alive_.size(); ++i) {
            if (alive_[i]) {
                level.originalToLocal[i] = localIdx;
                level.localToOriginal.push_back(i);
                level.positions.push_back(positions_[i]);
                localIdx++;
            }
        }

        // Remap faces
        for (std::size_t fi = 0; fi < faces_.size(); ++fi) {
            if (!faceAlive_[fi]) {
                continue;
            }
            std::array<std::size_t, 3> localTri;
            for (int j = 0; j < 3; ++j) {
                localTri[j] = level.originalToLocal.at(faces_[fi][j]);
            }
            level.faces.push_back(localTri);
        }

        return level;
    }

    /** Rebuild the edge list from alive faces and return it */
    auto rebuildAndGetEdges() -> const std::vector<std::pair<std::size_t, std::size_t>>&
    {
        buildEdges_();
        return edges_;
    }

private:
    void computeQuadrics_()
    {
        for (auto& q : quadrics_) {
            q = Quadric<T>();
        }

        for (std::size_t fi = 0; fi < faces_.size(); ++fi) {
            if (!faceAlive_[fi]) {
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

    void buildEdges_()
    {
        edges_.clear();
        std::unordered_set<std::size_t> seen;
        auto edgeKey = [this](std::size_t a, std::size_t b) -> std::size_t {
            auto n = positions_.size();
            return std::min(a, b) * n + std::max(a, b);
        };

        for (std::size_t fi = 0; fi < faces_.size(); ++fi) {
            if (!faceAlive_[fi]) {
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
        T u = T(1) - v - w;

        return {u, v, w};
    }

    std::vector<Vec<T, 3>> positions_;
    std::vector<bool> alive_;
    std::vector<bool> isBoundary_;
    std::vector<bool> isPinned_;
    std::vector<Quadric<T>> quadrics_;
    std::vector<std::array<std::size_t, 3>> faces_;
    std::vector<bool> faceAlive_;
    std::vector<std::vector<std::size_t>> vertFaces_;
    std::size_t numAliveVerts_{0};
    std::size_t numAliveFaces_{0};
    std::vector<std::pair<std::size_t, std::size_t>> edges_;
};

/**
 * @brief Build a mesh hierarchy by greedy QEM decimation
 *
 * Returns a vector of HierarchyLevel from finest to coarsest, plus
 * the collapse records needed for prolongation (ordered from finest to coarsest).
 */
template <typename T, class MeshPtr>
auto buildHierarchy(const MeshPtr& mesh, std::size_t pin0, std::size_t pin1, std::size_t levelRatio,
                    std::size_t minCoarseVerts)
    -> std::pair<std::vector<HierarchyLevel<T>>, std::vector<std::vector<CollapseRecord<T>>>>
{
    DecimationMesh<T> dmesh;
    dmesh.build(mesh, pin0, pin1);

    // Finest level snapshot
    std::vector<HierarchyLevel<T>> levels;
    levels.push_back(dmesh.snapshot());

    std::vector<std::vector<CollapseRecord<T>>> collapsesByLevel;

    auto targetVerts = dmesh.numAliveVerts();
    if (targetVerts <= minCoarseVerts) {
        // Mesh is already small enough — single level
        return {levels, collapsesByLevel};
    }

    while (targetVerts > minCoarseVerts) {
        auto nextTarget = std::max(targetVerts / levelRatio, minCoarseVerts);
        std::vector<CollapseRecord<T>> levelCollapses;

        // Build priority queue of edge collapses
        using CostEdge = std::pair<T, std::pair<std::size_t, std::size_t>>;
        std::priority_queue<CostEdge, std::vector<CostEdge>, std::greater<CostEdge>> pq;

        const auto& edges = dmesh.rebuildAndGetEdges();
        for (auto& [a, b] : edges) {
            // Try collapsing the collapsible vertex towards the other
            if (dmesh.isCollapsible(a)) {
                pq.push({dmesh.collapseCost(a, b), {a, b}});
            }
            if (dmesh.isCollapsible(b)) {
                pq.push({dmesh.collapseCost(b, a), {b, a}});
            }
        }

        while (dmesh.numAliveVerts() > nextTarget && !pq.empty()) {
            auto [cost, edge] = pq.top();
            pq.pop();

            auto [vRemove, vKeep] = edge;
            auto record = dmesh.tryCollapse(vRemove, vKeep);
            if (!record) {
                continue;
            }

            levelCollapses.push_back(*record);

            // Add new edges involving vKeep to the priority queue
            auto nbrs = dmesh.vertexNeighbors(vKeep);
            for (auto nb : nbrs) {
                if (dmesh.isCollapsible(nb)) {
                    pq.push({dmesh.collapseCost(nb, vKeep), {nb, vKeep}});
                }
                if (dmesh.isCollapsible(vKeep)) {
                    pq.push({dmesh.collapseCost(vKeep, nb), {vKeep, nb}});
                }
            }
        }

        if (levelCollapses.empty()) {
            break;  // No more valid collapses possible
        }

        collapsesByLevel.push_back(std::move(levelCollapses));
        levels.push_back(dmesh.snapshot());
        targetVerts = dmesh.numAliveVerts();
    }

    return {levels, collapsesByLevel};
}

/**
 * @brief Build a HalfEdgeMesh from a hierarchy level
 */
template <typename T>
auto buildLevelMesh(const HierarchyLevel<T>& level) -> typename HalfEdgeMesh<T>::Pointer
{
    auto mesh = HalfEdgeMesh<T>::New();

    for (const auto& pos : level.positions) {
        mesh->insert_vertex(pos[0], pos[1], pos[2]);
    }

    std::vector<std::vector<std::size_t>> faceVec;
    faceVec.reserve(level.faces.size());
    for (const auto& tri : level.faces) {
        faceVec.push_back({tri[0], tri[1], tri[2]});
    }
    mesh->insert_faces(faceVec);

    return mesh;
}

/**
 * @brief Prolongate UV coordinates from a coarser level to a finer level
 *
 * Surviving vertices get their UVs directly; removed vertices get UVs
 * via barycentric interpolation in their containing post-collapse triangle.
 *
 * @param coarseUVs UV coordinates indexed by original vertex index
 * @param collapses Collapse records for this level transition (finest-to-coarsest order)
 * @return UV map indexed by original vertex index (includes all finer-level vertices)
 */
template <typename T>
auto prolongateUVs(const std::unordered_map<std::size_t, std::array<T, 2>>& coarseUVs,
                   const std::vector<CollapseRecord<T>>& collapses)
    -> std::unordered_map<std::size_t, std::array<T, 2>>
{
    // Start with all coarse-level UVs
    auto fineUVs = coarseUVs;

    // Undo collapses in reverse order (coarsest collapse first was last applied)
    for (auto it = collapses.rbegin(); it != collapses.rend(); ++it) {
        auto& rec = *it;
        auto& tri = rec.containingTri;

        // All three containing-tri vertices should have UVs by now
        auto uv0 = fineUVs.at(tri[0]);
        auto uv1 = fineUVs.at(tri[1]);
        auto uv2 = fineUVs.at(tri[2]);

        std::array<T, 2> newUV;
        newUV[0] = rec.bary[0] * uv0[0] + rec.bary[1] * uv1[0] + rec.bary[2] * uv2[0];
        newUV[1] = rec.bary[0] * uv0[1] + rec.bary[1] * uv1[1] + rec.bary[2] * uv2[1];
        fineUVs[rec.vRemoved] = newUV;
    }

    return fineUVs;
}

/**
 * @brief Solve the LSCM system at one hierarchy level
 *
 * Builds the Lévy et al. Eq. 10 LSCM system on the given level mesh with the
 * given pin vertices. If an initial guess is provided, uses solveWithGuess.
 *
 * @return UV coordinates indexed by original vertex index
 */
template <typename T, class SolverType>
auto solveLSCMLevel(const typename HalfEdgeMesh<T>::Pointer& levelMesh,
                    const detail::hlscm::HierarchyLevel<T>& level, std::size_t origPin0,
                    std::size_t origPin1,
                    const std::unordered_map<std::size_t, std::array<T, 2>>* initialGuess)
    -> std::unordered_map<std::size_t, std::array<T, 2>>
{
    using Triplet = Eigen::Triplet<T>;
    using SparseMatrix = Eigen::SparseMatrix<T>;
    using DenseMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

    auto numFaces = levelMesh->num_faces();
    auto numVerts = levelMesh->num_vertices();

    // Map original pin indices to level-local indices
    auto localPin0 = level.originalToLocal.at(origPin0);
    auto localPin1 = level.originalToLocal.at(origPin1);
    auto p0 = levelMesh->vertex(localPin0);
    auto p1 = levelMesh->vertex(localPin1);

    // Place pins on UV axes (same logic as AngleBasedLSCM)
    auto pinVec = p1->pos - p0->pos;
    auto dist = norm(pinVec);
    pinVec /= dist;
    p0->pos = {T(0), T(0), T(0)};
    auto maxElem = std::max_element(pinVec.begin(), pinVec.end());
    auto maxAxis = std::distance(pinVec.begin(), maxElem);
    dist = std::copysign(dist, *maxElem);
    if (maxAxis == 0) {
        p1->pos = {dist, T(0), T(0)};
    } else {
        p1->pos = {T(0), dist, T(0)};
    }

    auto numFixed = std::size_t(2);
    auto numFree = numVerts - numFixed;

    // Build free vertex index table
    std::map<std::size_t, std::size_t> freeIdxTable;
    for (const auto& v : levelMesh->vertices()) {
        if (v == p0 || v == p1) {
            continue;
        }
        auto newIdx = freeIdxTable.size();
        freeIdxTable[v->idx] = newIdx;
    }

    // Setup pinned bFixed
    std::vector<Triplet> tripletsB;
    tripletsB.emplace_back(0, 0, p0->pos[0]);
    tripletsB.emplace_back(1, 0, p0->pos[1]);
    tripletsB.emplace_back(2, 0, p1->pos[0]);
    tripletsB.emplace_back(3, 0, p1->pos[1]);
    SparseMatrix bFixed(2 * numFixed, 1);
    bFixed.reserve(tripletsB.size());
    bFixed.setFromTriplets(tripletsB.begin(), tripletsB.end());

    // Build LSCM system matrices
    std::vector<Triplet> tripletsA;
    tripletsB.clear();

    auto addContrib = [&](std::size_t row, const auto& e, T c, T s) {
        if (e->vertex == p0) {
            tripletsB.emplace_back(row, 0, c);
            tripletsB.emplace_back(row, 1, -s);
            tripletsB.emplace_back(row + 1, 0, s);
            tripletsB.emplace_back(row + 1, 1, c);
        } else if (e->vertex == p1) {
            tripletsB.emplace_back(row, 2, c);
            tripletsB.emplace_back(row, 3, -s);
            tripletsB.emplace_back(row + 1, 2, s);
            tripletsB.emplace_back(row + 1, 3, c);
        } else {
            auto freeIdx = freeIdxTable.at(e->vertex->idx);
            tripletsA.emplace_back(row, 2 * freeIdx, c);
            tripletsA.emplace_back(row, 2 * freeIdx + 1, -s);
            tripletsA.emplace_back(row + 1, 2 * freeIdx, s);
            tripletsA.emplace_back(row + 1, 2 * freeIdx + 1, c);
        }
    };

    for (const auto& f : levelMesh->faces()) {
        auto e0 = f->head;
        auto e1 = e0->next;
        auto e2 = e1->next;
        auto sin0 = std::sin(e0->alpha);
        auto sin1 = std::sin(e1->alpha);
        auto sin2 = std::sin(e2->alpha);

        std::vector<T> sins{sin0, sin1, sin2};
        auto sinMaxElem = std::max_element(sins.begin(), sins.end());
        auto sinMaxIdx = std::distance(sins.begin(), sinMaxElem);

        if (sinMaxIdx == 0) {
            auto temp = e0;
            e0 = e1;
            e1 = e2;
            e2 = temp;
            sin0 = sins[1];
            sin1 = sins[2];
            sin2 = sins[0];
        } else if (sinMaxIdx == 1) {
            auto temp = e2;
            e2 = e1;
            e1 = e0;
            e0 = temp;
            sin0 = sins[2];
            sin1 = sins[0];
            sin2 = sins[1];
        }

        auto ratio = (sin2 == T(0)) ? T(1) : sin1 / sin2;
        auto cosine = std::cos(e0->alpha) * ratio;
        auto sine = sin0 * ratio;

        auto row = 2 * f->idx;
        addContrib(row, e0, cosine - T(1), sine);
        addContrib(row, e1, -cosine, -sine);
        addContrib(row, e2, T(1), T(0));
    }

    SparseMatrix A(2 * numFaces, 2 * numFree);
    A.reserve(tripletsA.size());
    A.setFromTriplets(tripletsA.begin(), tripletsA.end());

    SparseMatrix bFree(2 * numFaces, 2 * numFixed);
    bFree.reserve(tripletsB.size());
    bFree.setFromTriplets(tripletsB.begin(), tripletsB.end());

    SparseMatrix b = bFree * bFixed * T(-1);

    // Solve
    DenseMatrix x;
    if constexpr (detail::is_instance_of_v<SolverType, Eigen::LeastSquaresConjugateGradient>) {
        // LSCG solves the rectangular system A x = b directly
        if (initialGuess && !initialGuess->empty()) {
            // Build initial guess vector from prolongated UVs
            DenseMatrix x0(2 * numFree, 1);
            for (const auto& v : levelMesh->vertices()) {
                if (v == p0 || v == p1) {
                    continue;
                }
                auto freeIdx = freeIdxTable.at(v->idx);
                auto origIdx = level.localToOriginal[v->idx];
                auto guessIt = initialGuess->find(origIdx);
                if (guessIt != initialGuess->end()) {
                    x0(2 * freeIdx, 0) = guessIt->second[0];
                    x0(2 * freeIdx + 1, 0) = guessIt->second[1];
                } else {
                    x0(2 * freeIdx, 0) = T(0);
                    x0(2 * freeIdx + 1, 0) = T(0);
                }
            }
            DenseMatrix bDense = b;
            SolverType solver(A);
            x = solver.solveWithGuess(bDense, x0);
            if (solver.info() != Eigen::ComputationInfo::Success) {
                throw SolverException("HLSCM: LSCG solve failed at hierarchy level");
            }
        } else {
            SolverType solver(A);
            DenseMatrix bDense = b;
            x = solver.solve(bDense);
            if (solver.info() != Eigen::ComputationInfo::Success) {
                throw SolverException("HLSCM: LSCG solve failed at hierarchy level");
            }
        }
    } else if constexpr (std::is_base_of_v<Eigen::IterativeSolverBase<SolverType>, SolverType>) {
        // Other iterative solvers (e.g. ConjugateGradient) require a square SPD matrix;
        // use normal equations AtA x = Atb
        SparseMatrix AtA = A.transpose() * A;
        AtA.makeCompressed();
        SparseMatrix Atb = A.transpose() * b;
        if (initialGuess && !initialGuess->empty()) {
            // Build initial guess vector from prolongated UVs
            DenseMatrix x0(2 * numFree, 1);
            for (const auto& v : levelMesh->vertices()) {
                if (v == p0 || v == p1) {
                    continue;
                }
                auto freeIdx = freeIdxTable.at(v->idx);
                auto origIdx = level.localToOriginal[v->idx];
                auto guessIt = initialGuess->find(origIdx);
                if (guessIt != initialGuess->end()) {
                    x0(2 * freeIdx, 0) = guessIt->second[0];
                    x0(2 * freeIdx + 1, 0) = guessIt->second[1];
                } else {
                    x0(2 * freeIdx, 0) = T(0);
                    x0(2 * freeIdx + 1, 0) = T(0);
                }
            }
            DenseMatrix AtbDense = Atb;
            SolverType solver(AtA);
            x = solver.solveWithGuess(AtbDense, x0);
            if (solver.info() != Eigen::ComputationInfo::Success) {
                throw SolverException("HLSCM: iterative solve failed at hierarchy level");
            }
        } else {
            DenseMatrix AtbDense = Atb;
            SolverType solver(AtA);
            x = solver.solve(AtbDense);
            if (solver.info() != Eigen::ComputationInfo::Success) {
                throw SolverException("HLSCM: iterative solve failed at hierarchy level");
            }
        }
    } else {
        // Direct solver: no initial guess support; ignore initialGuess
        x = detail::SolveLeastSquares<SparseMatrix, DenseMatrix, SolverType>(A, b);
    }

    // Build output UV map (original vertex indices → UV)
    std::unordered_map<std::size_t, std::array<T, 2>> uvs;
    uvs[level.localToOriginal[p0->idx]] = {p0->pos[0], p0->pos[1]};
    uvs[level.localToOriginal[p1->idx]] = {p1->pos[0], p1->pos[1]};
    for (const auto& v : levelMesh->vertices()) {
        if (v == p0 || v == p1) {
            continue;
        }
        auto freeIdx = 2 * freeIdxTable.at(v->idx);
        auto origIdx = level.localToOriginal[v->idx];
        uvs[origIdx] = {x(freeIdx, 0), x(freeIdx + 1, 0)};
    }
    return uvs;
}

}  // namespace hlscm
}  // namespace detail

/**
 * @brief Compute parameterized mesh using Hierarchical LSCM
 *
 * Implements the HLSCM algorithm from Ray & Lévy, "Hierarchical Least Squares
 * Conformal Map" (2003). Uses cascadic multigrid to accelerate LSCM
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
 *         initial guess. Defaults to LeastSquaresConjugateGradient, which
 *         operates directly on the rectangular system and yields the best
 *         convergence rate when combined with the hierarchical warm-start.
 *         ConjugateGradient can be used for flat LSCM (no hierarchy) where
 *         it is faster, but provides no benefit inside HLSCM.
 */
template <typename T, class MeshType = HalfEdgeMesh<T>,
          class Solver = Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<T>>,
          std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
class HierarchicalLSCM
{
public:
    /** @brief Mesh type alias */
    using Mesh = MeshType;

    /** @brief Set the pinned vertex indices used by compute() */
    void setPinnedVertices(std::size_t pin0Idx, std::size_t pin1Idx)
    {
        pinnedVertices_ = {pin0Idx, pin1Idx};
    }

    /** @brief Set the vertex ratio between consecutive hierarchy levels */
    void setLevelRatio(std::size_t ratio) { levelRatio_ = ratio; }

    /** @brief Set the minimum vertex count at the coarsest level */
    void setMinCoarseVertices(std::size_t count) { minCoarseVertices_ = count; }

    /** @copydoc HierarchicalLSCM::Compute() */
    void compute(typename Mesh::Pointer& mesh) const
    {
        std::size_t p0, p1;
        if (pinnedVertices_) {
            p0 = pinnedVertices_->first;
            p1 = pinnedVertices_->second;
        } else {
            AutoSelectPins(mesh, p0, p1);
        }
        ComputeImpl(mesh, p0, p1, levelRatio_, minCoarseVertices_);
    }

    /**
     * @brief Compute with automatic pin selection
     *
     * Selects pins identically to AngleBasedLSCM::Compute().
     */
    static void Compute(typename Mesh::Pointer& mesh)
    {
        std::size_t p0, p1;
        AutoSelectPins(mesh, p0, p1);
        ComputeImpl(mesh, p0, p1);
    }

    /**
     * @brief Compute with explicit pinned vertex indices
     */
    static void Compute(typename Mesh::Pointer& mesh, std::size_t pin0Idx, std::size_t pin1Idx)
    {
        ComputeImpl(mesh, pin0Idx, pin1Idx);
    }

private:
    /** Select two pinned boundary vertices (same logic as AngleBasedLSCM) */
    static void AutoSelectPins(const typename Mesh::Pointer& mesh, std::size_t& p0, std::size_t& p1)
    {
        auto v0 = mesh->vertices_boundary().front();
        auto e = v0->edge;
        do {
            if (e->pair->is_boundary()) {
                break;
            }
            e = e->pair->next;
        } while (e != v0->edge);
        if (e == v0->edge && !e->pair->is_boundary()) {
            throw MeshException("Pinned vertex not on boundary");
        }
        p0 = v0->idx;
        p1 = e->next->vertex->idx;
    }

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

    static void ComputeImpl(typename Mesh::Pointer& mesh, std::size_t pin0Idx, std::size_t pin1Idx,
                            std::size_t levelRatio = 10, std::size_t minCoarseVerts = 100)
    {
        // Build mesh hierarchy
        auto [levels, collapsesByLevel] =
            detail::hlscm::buildHierarchy<T>(mesh, pin0Idx, pin1Idx, levelRatio, minCoarseVerts);

        if (levels.size() <= 1) {
            // Mesh too small for hierarchy — single-level LSCM solve
            // Use AngleBasedLSCM for exact equivalence on small meshes
            AngleBasedLSCM<T, MeshType, Solver>::Compute(mesh, pin0Idx, pin1Idx);
            return;
        }

        // Solve coarsest level (last in the array)
        auto coarsestIdx = levels.size() - 1;
        auto coarseMesh = detail::hlscm::buildLevelMesh<T>(levels[coarsestIdx]);
        ComputeMeshAngles(coarseMesh);
        auto uvs = detail::hlscm::solveLSCMLevel<T, Solver>(coarseMesh, levels[coarsestIdx],
                                                            pin0Idx, pin1Idx, nullptr);

        // Prolongate and refine at each finer level
        for (std::size_t k = coarsestIdx; k-- > 0;) {
            // Prolongate UVs from level k+1 to level k
            uvs = detail::hlscm::prolongateUVs<T>(uvs, collapsesByLevel[k]);

            // Build level mesh
            auto levelMesh = detail::hlscm::buildLevelMesh<T>(levels[k]);

            if (k == 0) {
                // Finest level: use original mesh angles (may be ABF-optimized)
                CopyAnglesFromOriginal(mesh, levelMesh);
            } else {
                // Coarser levels: compute angles from 3D geometry
                ComputeMeshAngles(levelMesh);
            }

            // Solve with initial guess
            uvs = detail::hlscm::solveLSCMLevel<T, Solver>(levelMesh, levels[k], pin0Idx, pin1Idx,
                                                           &uvs);
        }

        // Transfer final UVs back to input mesh
        for (const auto& v : mesh->vertices()) {
            auto it = uvs.find(v->idx);
            if (it != uvs.end()) {
                v->pos[0] = it->second[0];
                v->pos[1] = it->second[1];
                v->pos[2] = T(0);
            }
        }
    }

    /** Optional explicit pin pair */
    std::optional<std::pair<std::size_t, std::size_t>> pinnedVertices_;
    std::size_t levelRatio_{10};
    std::size_t minCoarseVertices_{100};
};

}  // namespace OpenABF
