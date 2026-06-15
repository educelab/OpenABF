#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <Eigen/SparseCore>

#include "OpenABF/Exceptions.hpp"
#include "OpenABF/Math.hpp"
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
 * @brief Validate a user-supplied PinMap against `mesh`.
 *
 * Throws `std::invalid_argument` if the PinMap has fewer than two pins, a
 * duplicate vertex index, or an out-of-range vertex index. Shared by
 * `AngleBasedLSCM` and `HierarchicalLSCM` so behavior matches at the public
 * boundary.
 */
template <typename T, class MeshType>
void ValidatePins(const typename MeshType::Pointer& mesh, const PinMap<T>& pins)
{
    if (pins.size() < 2) {
        throw std::invalid_argument("LSCM: PinMap requires at least 2 pins");
    }
    const auto numVerts = mesh->num_vertices();
    std::unordered_set<std::size_t> seen;
    seen.reserve(pins.size());
    for (const auto& [vIdx, uv] : pins) {
        if (vIdx >= numVerts) {
            throw std::invalid_argument("LSCM: PinMap vertex index out of range");
        }
        if (!seen.insert(vIdx).second) {
            throw std::invalid_argument("LSCM: PinMap has duplicate vertex index");
        }
    }
}

/**
 * @brief Build a 2-entry PinMap from explicit vertex indices using the LSCM
 * axis-snap convention.
 *
 * pin0 lands at `{0, 0}`; pin1 lands at signed distance `|p1 - p0|` on
 * whichever world-axis the `(p1 - p0)` vector has the largest magnitude.
 * Used by the auto-pin path and by the deprecated 2-pin shims.
 *
 * Throws `std::invalid_argument` if `p0Idx == p1Idx` (the resulting
 * zero-length axis-snap would divide by zero and produce NaN UVs) or if
 * either index is out of range.
 */
template <typename T, class MeshType>
auto AutoPlacePair(const typename MeshType::Pointer& mesh, std::size_t p0Idx,
                   std::size_t p1Idx) -> PinMap<T>
{
    const auto numVerts = mesh->num_vertices();
    if (p0Idx >= numVerts || p1Idx >= numVerts) {
        throw std::invalid_argument("LSCM: pin vertex index out of range");
    }
    if (p0Idx == p1Idx) {
        throw std::invalid_argument("LSCM: pin pair must be two distinct vertices");
    }
    auto p0 = mesh->vertex(p0Idx);
    auto p1 = mesh->vertex(p1Idx);
    auto pinVec = p1->pos - p0->pos;
    auto dist = norm(pinVec);
    pinVec /= dist;
    auto maxElem = std::max_element(pinVec.begin(), pinVec.end());
    auto maxAxis = std::distance(pinVec.begin(), maxElem);
    dist = std::copysign(dist, *maxElem);
    Vec<T, 2> uv0{T(0), T(0)};
    Vec<T, 2> uv1 = (maxAxis == 0) ? Vec<T, 2>{dist, T(0)} : Vec<T, 2>{T(0), dist};
    return PinMap<T>{{p0Idx, uv0}, {p1Idx, uv1}};
}

/**
 * @brief Auto-select two boundary pins and place them via the LSCM
 * axis-snap convention.
 *
 * Picks the first boundary vertex returned by `mesh->vertices_boundary()` as
 * pin0 and walks the boundary to find an adjacent boundary vertex as pin1.
 * UVs follow the axis-snap convention applied by `AutoPlacePair`.
 *
 * @throws MeshException if the mesh has no boundary vertices, or if no
 *         boundary-adjacent neighbor is found for the first boundary vertex.
 */
template <typename T, class MeshType>
auto AutoSelectPins(const typename MeshType::Pointer& mesh) -> PinMap<T>
{
    auto boundary = mesh->vertices_boundary();
    if (boundary.empty()) {
        throw MeshException("LSCM: mesh has no boundary vertices");
    }
    auto p0 = boundary.front();
    auto e = p0->edge;
    do {
        if (e->pair->is_boundary()) {
            break;
        }
        e = e->pair->next;
    } while (e != p0->edge);
    if (e == p0->edge && !e->pair->is_boundary()) {
        throw MeshException("LSCM: pinned vertex not on boundary");
    }
    auto p1 = e->next->vertex;
    return AutoPlacePair<T, MeshType>(mesh, p0->idx, p1->idx);
}

/**
 * @brief Outputs of `BuildSystem`: the LSCM least-squares system for a mesh
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
 * Pure with respect to the mesh — only reads `alpha` and connectivity. Pin
 * UVs are taken verbatim from the PinMap into `bFixed`; the mesh's vertex
 * positions are NOT mutated. Callers that need pin UVs reflected on the mesh
 * (e.g., `AngleBasedLSCM::ComputeImpl`'s output writeback) must do that
 * themselves.
 *
 * The function does NOT auto-place pins — UVs are taken verbatim from the
 * PinMap. Callers that want the LSCM axis-snap convention (origin +
 * dominant-axis placement for an auto-selected pair) compute those UVs
 * themselves via `AutoPlacePair` before calling this helper.
 *
 * Assembly follows Lévy et al. 2002 Eq. 10 using the per-edge `alpha` angles
 * already stored on the mesh.
 */
template <typename T, class MeshType>
auto BuildSystem(const typename MeshType::Pointer& mesh, const PinMap<T>& pins) -> SystemParts<T>
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

    // Populate bFixed from PinMap UVs (verbatim).
    std::vector<Triplet> tripletsB;
    tripletsB.reserve(2 * numFixed);
    for (std::size_t s = 0; s < numFixed; ++s) {
        const auto& uv = pins[s].second;
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
