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
    /** Index of the third vertex of the containing triangle (for barycentric interp) */
    std::size_t vThird;
    /** Barycentric coordinates (w.r.t. vKept, vThird, vRemoved) of vRemoved */
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

        // Check angular defect: ensure no face would flip after collapse
        for (auto fi : removeFaces) {
            // This face contains vRemove but not vKeep.
            // After collapse, vRemove is replaced by vKeep.
            std::array<std::size_t, 3> newTri = faces_[fi];
            for (auto& vi : newTri) {
                if (vi == vRemove) {
                    vi = vKeep;
                }
            }
            // Check that the triangle normal doesn't flip
            auto& p0 = positions_[newTri[0]];
            auto& p1 = positions_[newTri[1]];
            auto& p2 = positions_[newTri[2]];
            auto e1 = p1 - p0;
            auto e2 = p2 - p0;
            auto newNormal = cross(e1, e2);
            auto newArea = norm(newNormal);
            if (newArea < std::numeric_limits<T>::epsilon() * T(100)) {
                return std::nullopt;  // Degenerate triangle
            }

            // Check against old normal
            auto& op0 = positions_[faces_[fi][0]];
            auto& op1 = positions_[faces_[fi][1]];
            auto& op2 = positions_[faces_[fi][2]];
            auto oe1 = op1 - op0;
            auto oe2 = op2 - op0;
            auto oldNormal = cross(oe1, oe2);
            if (dot(oldNormal, newNormal) < T(0)) {
                return std::nullopt;  // Normal flip
            }
        }

        // Build collapse record with barycentric coordinates
        // vRemoved's position in terms of the nearest surviving triangle
        // For prolongation, we store bary coords of vRemove w.r.t. a
        // triangle containing vKeep. We use one of the updated removeFaces.
        CollapseRecord<T> record;
        record.vRemoved = vRemove;
        record.vKept = vKeep;

        if (!removeFaces.empty()) {
            // Use the first surviving face that will contain vKeep after collapse
            auto fi = removeFaces[0];
            std::size_t third = 0;
            for (auto vi : faces_[fi]) {
                if (vi != vRemove && vi != vKeep) {
                    third = vi;
                    break;
                }
            }
            record.vThird = third;
            record.bary = computeBarycentric_(positions_[vRemove], positions_[vKeep],
                                              positions_[third], positions_[vRemove]);
        } else {
            // Edge case: all faces are shared (degree-2 vertex) — shouldn't happen
            // for interior vertex with 2 shared faces, but handle gracefully
            record.vThird = vKeep;
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

    /** Get all edges as directed pairs (v0 < v1 to avoid duplicates) */
    [[nodiscard]] auto allEdges() const -> std::vector<std::pair<std::size_t, std::size_t>>
    {
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

        auto edges = dmesh.allEdges();
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
 * @tparam Solver An Eigen iterative solver supporting solveWithGuess
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
        if (pinnedVertices_) {
            Compute(mesh, pinnedVertices_->first, pinnedVertices_->second);
        } else {
            Compute(mesh);
        }
    }

    /**
     * @brief Compute with automatic pin selection
     *
     * Selects pins identically to AngleBasedLSCM::Compute().
     */
    static void Compute(typename Mesh::Pointer& mesh)
    {
        // Delegate to AngleBasedLSCM for now (stub)
        AngleBasedLSCM<T, MeshType, Solver>::Compute(mesh);
    }

    /**
     * @brief Compute with explicit pinned vertex indices
     */
    static void Compute(typename Mesh::Pointer& mesh, std::size_t pin0Idx, std::size_t pin1Idx)
    {
        // Delegate to AngleBasedLSCM for now (stub)
        AngleBasedLSCM<T, MeshType, Solver>::Compute(mesh, pin0Idx, pin1Idx);
    }

private:
    /** Optional explicit pin pair */
    std::optional<std::pair<std::size_t, std::size_t>> pinnedVertices_;
    std::size_t levelRatio_{10};
    std::size_t minCoarseVertices_{100};
};

}  // namespace OpenABF
