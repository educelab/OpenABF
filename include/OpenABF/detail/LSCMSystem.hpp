#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <unordered_map>
#include <vector>

#include <Eigen/SparseCore>

namespace OpenABF::detail::lscm
{

/**
 * @internal
 * @brief Outputs of `buildSystem`: the LSCM least-squares system for a mesh
 *        with two pinned vertices.
 *
 * @tparam T Floating-point type
 *
 * Layout: `A` is `(2·numFaces) × (2·numFree)`, `b` is `(2·numFaces) × 1`,
 * `freeIdxTable` maps `vertex->idx` to a row-pair index in `A`/`x` (so a free
 * vertex `v` occupies rows `2*freeIdxTable[v->idx]` and
 * `2*freeIdxTable[v->idx]+1`).
 */
template <typename T>
struct SystemParts {
    /** Coefficient matrix of the LSCM least-squares system. Shape `(2·numFaces) × (2·numFree)`. */
    Eigen::SparseMatrix<T> A;
    /** Right-hand side vector. Shape `(2·numFaces) × 1`. Contains the pin contributions. */
    Eigen::SparseMatrix<T> b;
    /** Maps mesh vertex `idx` to a row-pair slot in `A`/`x`. Excludes pin vertices. */
    std::unordered_map<std::size_t, std::size_t> freeIdxTable;
};

/**
 * @internal
 * @brief Build the LSCM sparse system for a mesh with two pinned vertices.
 *
 * Mutates the mesh: places `p0` at the UV origin and `p1` on whichever XY
 * axis its displacement from `p0` has the largest magnitude — same pin
 * placement convention used by `AngleBasedLSCM::ComputeImpl` and
 * `HierarchicalLSCM::solveLSCMLevel`.
 *
 * Assembly follows Lévy et al. 2002 Eq. 10 using the per-edge `alpha` angles
 * already stored on the mesh.
 */
template <typename T, class MeshType>
auto buildSystem(const typename MeshType::Pointer& mesh, const typename MeshType::VertPtr& p0,
                 const typename MeshType::VertPtr& p1) -> SystemParts<T>
{
    using Triplet = Eigen::Triplet<T>;
    using SparseMatrix = Eigen::SparseMatrix<T>;

    // Map selected edge to closest XY axis. Use sign to select direction.
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

    const auto numFaces = mesh->num_faces();
    const auto numVerts = mesh->num_vertices();
    constexpr std::size_t numFixed = 2;
    const auto numFree = numVerts - numFixed;

    // Permutation for free vertices: maps mesh vertex idx → row-pair slot in A.
    std::unordered_map<std::size_t, std::size_t> freeIdxTable;
    freeIdxTable.reserve(numFree);
    for (const auto& v : mesh->vertices()) {
        if (v == p0 or v == p1) {
            continue;
        }
        auto newIdx = freeIdxTable.size();
        freeIdxTable[v->idx] = newIdx;
    }

    // Setup pinned bFixed.
    std::vector<Triplet> tripletsB;
    tripletsB.emplace_back(0, 0, p0->pos[0]);
    tripletsB.emplace_back(1, 0, p0->pos[1]);
    tripletsB.emplace_back(2, 0, p1->pos[0]);
    tripletsB.emplace_back(3, 0, p1->pos[1]);
    SparseMatrix bFixed(2 * numFixed, 1);
    bFixed.reserve(tripletsB.size());
    bFixed.setFromTriplets(tripletsB.begin(), tripletsB.end());

    // Setup variables matrix. Only solving for free vertices, so pins go in
    // a special matrix.
    std::vector<Triplet> tripletsA;
    tripletsB.clear();

    // Per-vertex contribution helper (Lévy et al. 2002, Eq. 10).
    // Each vertex contributes a 2×2 conformal block [c, -s; s, c] at its
    // column. Fixed pins (p0, p1) go into tripletsB; free vertices into
    // tripletsA.
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

    for (const auto& f : mesh->faces()) {
        auto e0 = f->head;
        auto e1 = e0->next;
        auto e2 = e1->next;
        auto sin0 = std::sin(e0->alpha);
        auto sin1 = std::sin(e1->alpha);
        auto sin2 = std::sin(e2->alpha);

        // Find the max sin idx and rotate the edge order so last angle is largest.
        std::array<T, 3> sins{sin0, sin1, sin2};
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

    return SystemParts<T>{std::move(A), std::move(b), std::move(freeIdxTable)};
}

}  // namespace OpenABF::detail::lscm
