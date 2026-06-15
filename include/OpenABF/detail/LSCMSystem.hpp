#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <unordered_map>
#include <utility>
#include <vector>

#include <Eigen/SparseCore>

#include "OpenABF/Vec.hpp"

namespace OpenABF::detail::lscm
{

/**
 * @brief Per-pin entry: (mesh vertex index, target UV).
 *
 * Shared by `AngleBasedLSCM::PinMap` and `HierarchicalLSCM::PinMap`.
 */
template <typename T>
using PinMap = std::vector<std::pair<std::size_t, OpenABF::Vec<T, 2>>>;

/**
 * @brief Outputs of `buildSystem`: the LSCM least-squares system for a mesh
 *        with N pinned vertices (N ≥ 2).
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
 * @brief Build the LSCM sparse system for a mesh with N pinned vertices.
 *
 * Mutates the mesh: each pinned vertex's `pos` is overwritten with `{uv[0],
 * uv[1], 0}` (the caller's chosen UV). The function does NOT auto-place pins
 * — UVs are taken verbatim from the PinMap. Callers that want the LSCM
 * axis-snap convention (origin + dominant-axis placement for an auto-selected
 * pair) compute those UVs themselves before calling this helper.
 *
 * Assembly follows Lévy et al. 2002 Eq. 10 using the per-edge `alpha` angles
 * already stored on the mesh.
 */
template <typename T, class MeshType>
auto buildSystem(const typename MeshType::Pointer& mesh, const PinMap<T>& pins) -> SystemParts<T>
{
    using Triplet = Eigen::Triplet<T>;
    using SparseMatrix = Eigen::SparseMatrix<T>;

    const auto numFaces = mesh->num_faces();
    const auto numVerts = mesh->num_vertices();
    const auto numFixed = pins.size();
    const auto numFree = numVerts - numFixed;

    // Pin index → slot (0-based position within the PinMap).
    std::unordered_map<std::size_t, std::size_t> pinSlot;
    pinSlot.reserve(numFixed);
    for (std::size_t s = 0; s < numFixed; ++s) {
        pinSlot.emplace(pins[s].first, s);
    }

    // Write pin UVs into the mesh and into bFixed.
    std::vector<Triplet> tripletsB;
    tripletsB.reserve(2 * numFixed);
    for (std::size_t s = 0; s < numFixed; ++s) {
        const auto& [vIdx, uv] = pins[s];
        auto v = mesh->vertex(vIdx);
        v->pos = {uv[0], uv[1], T(0)};
        tripletsB.emplace_back(2 * s, 0, uv[0]);
        tripletsB.emplace_back(2 * s + 1, 0, uv[1]);
    }
    SparseMatrix bFixed(2 * numFixed, 1);
    bFixed.reserve(tripletsB.size());
    bFixed.setFromTriplets(tripletsB.begin(), tripletsB.end());

    // Permutation for free vertices: maps mesh vertex idx → row-pair slot in A.
    std::unordered_map<std::size_t, std::size_t> freeIdxTable;
    freeIdxTable.reserve(numFree);
    for (const auto& v : mesh->vertices()) {
        if (pinSlot.count(v->idx)) {
            continue;
        }
        auto newIdx = freeIdxTable.size();
        freeIdxTable[v->idx] = newIdx;
    }

    // Setup variables matrix. Only solving for free vertices, so pins go in
    // a special matrix.
    std::vector<Triplet> tripletsA;
    tripletsB.clear();

    // Per-vertex contribution helper (Lévy et al. 2002, Eq. 10).
    // Each vertex contributes a 2×2 conformal block [c, -s; s, c] at its
    // column. Pin vertices go into tripletsB at columns 2*slot and 2*slot+1;
    // free vertices into tripletsA.
    auto addContrib = [&](std::size_t row, const auto& e, T c, T s) {
        auto it = pinSlot.find(e->vertex->idx);
        if (it != pinSlot.end()) {
            auto col = 2 * it->second;
            tripletsB.emplace_back(row, col, c);
            tripletsB.emplace_back(row, col + 1, -s);
            tripletsB.emplace_back(row + 1, col, s);
            tripletsB.emplace_back(row + 1, col + 1, c);
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
