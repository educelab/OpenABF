#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <optional>
#include <type_traits>
#include <unordered_set>
#include <utility>

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>

#include "OpenABF/Exceptions.hpp"
#include "OpenABF/HalfEdgeMesh.hpp"
#include "OpenABF/Math.hpp"
#include "OpenABF/Vec.hpp"
#include "OpenABF/detail/LSCMSystem.hpp"

namespace OpenABF
{

namespace detail
{
/** Check if type is an instance of a template type: False */
template <class T, template <class...> class U>
constexpr bool is_instance_of_v = std::false_type{};

/** Check if type is an instance of a template type: True */
template <template <class...> class U, class... Vs>
constexpr bool is_instance_of_v<U<Vs...>, U> = std::true_type{};

/** Solve least squares using A'Ab  */
template <
    class SparseMatrix, class DenseMatrix, class Solver,
    std::enable_if_t<!is_instance_of_v<Solver, Eigen::LeastSquaresConjugateGradient>, bool> = false>
auto SolveLeastSquares(SparseMatrix A, SparseMatrix b) -> DenseMatrix
{
    // Setup AtA and solver
    SparseMatrix AtA = A.transpose() * A;
    AtA.makeCompressed();
    Solver solver;
    solver.compute(AtA);
    if (solver.info() != Eigen::ComputationInfo::Success) {
        throw SolverException("AB-LSCM: Failed to solve AtA");
    }

    // Setup Atb
    SparseMatrix Atb = A.transpose() * b;

    // Solve AtAx = AtAb
    DenseMatrix x = solver.solve(Atb);

    return x;
}

/** Solve least squares with LeastSquaresConjugateGradient */
template <
    class SparseMatrix, class DenseMatrix, class Solver,
    std::enable_if_t<is_instance_of_v<Solver, Eigen::LeastSquaresConjugateGradient>, bool> = true>
auto SolveLeastSquares(SparseMatrix A, SparseMatrix b) -> DenseMatrix
{
    // Solve
    Solver solver(A);
    DenseMatrix x = solver.solve(b);
    if (solver.info() != Eigen::ComputationInfo::Success) {
        throw SolverException("AB-LSCM: Failed to solve for b");
    }

    return x;
}

}  // namespace detail

/**
 * @brief Compute parameterized mesh using Angle-based LSCM
 *
 * Computes a least-squares conformal parameterization of a mesh. Unlike the
 * original LSCM algorithm, this class ignores the 3D vertex positions and
 * instead uses the angle associated with the mesh's edge trait
 * (MeshType::EdgeTraits::alpha) to calculate the initial per-triangle edge
 * lengths. Without previously modifying the angles of the provided mesh, this
 * class produces the same result as a vertex-based LSCM implementation.
 * However, by first processing the mesh with a parameterized angle optimizer,
 * such as ABFPlusPlus, the parameterization can be improved, sometimes
 * significantly.
 *
 * Implements the angle-based variant of "Least squares conformal maps for
 * automatic texture atlas generation" by Lévy _et al._ (2002)
 * \cite levy2002lscm.
 *
 * @tparam T Floating-point type
 * @tparam MeshType HalfEdgeMesh type which implements the default mesh traits
 * @tparam Solver A solver implementing the
 * [Eigen Sparse solver
 * concept](https://eigen.tuxfamily.org/dox-devel/group__TopicSparseSystems.html)
 * and templated on Eigen::SparseMatrix<T>. The default SparseLU is robust and
 * the fastest option we have benchmarked on the LSCM normal equations across
 * 50k–200k-face meshes (see `examples/src/BenchmarkFlattening.cpp`); for very
 * large meshes where SparseLU exhausts memory, switch to an iterative solver
 * via this template parameter.
 *
 * For iterative solving, prefer
 * `Eigen::ConjugateGradient<Eigen::SparseMatrix<T>, Eigen::Lower|Eigen::Upper>`
 * over the default `Lower`-only variant: the `Lower|Upper` template argument
 * enables Eigen's full-matrix SpMV code path, which is faster and — when
 * compiled with OpenMP — multi-threaded. Using only `Lower` (the Eigen
 * default) routes through `selfadjointView<Lower>`, which is a different
 * internal code path that is never OpenMP-parallelized regardless of
 * `Eigen::setNbThreads()`.
 *
 * `IncompleteCholesky` is also available as a preconditioner via
 * `Eigen::ConjugateGradient<..., Eigen::IncompleteCholesky<T>>`, but in our
 * benchmarks the per-iteration overhead of applying IC on the LSCM normal
 * equations dwarfs the iteration-count savings it provides over the default
 * `DiagonalPreconditioner` (Jacobi), and both are much slower than SparseLU
 * at the mesh sizes we target. Use IC only if you have profiled it favorably
 * against Diagonal for your specific mesh class.
 */
template <typename T, class MeshType = HalfEdgeMesh<T>,
          class Solver = Eigen::SparseLU<Eigen::SparseMatrix<T>, Eigen::COLAMDOrdering<int>>,
          std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
class AngleBasedLSCM
{
public:
    /** @brief Mesh type alias */
    using Mesh = MeshType;

    /**
     * @brief Per-pin entry: (mesh vertex index, target UV).
     *
     * A `PinMap` of size N ≥ 2 specifies an explicit LSCM pin set; each pinned
     * vertex's final UV equals its entry in the map.
     */
    using PinMap = detail::lscm::PinMap<T>;

    /**
     * @brief Set the explicit pin set used by `compute()`
     *
     * The PinMap must contain at least two unique, in-range vertex indices;
     * `compute()` rejects malformed inputs at solve time.
     */
    void setPins(PinMap pins)
    {
        pins_ = std::move(pins);
        legacyPinIndices_.reset();
    }

    /**
     * @brief Deprecated: set a pin pair by index using the LSCM axis-snap
     * convention at compute time.
     *
     * @deprecated Prefer `setPins(PinMap)`. This overload will be removed in
     * version 3.0.
     */
    [[deprecated("Use setPins(PinMap); will be removed in 3.0")]] void setPinnedVertices(
        std::size_t pin0Idx, std::size_t pin1Idx)
    {
        legacyPinIndices_ = {pin0Idx, pin1Idx};
        pins_.reset();
    }

    /** @copydoc AngleBasedLSCM::Compute() */
    void compute(typename Mesh::Pointer& mesh) const
    {
        if (pins_) {
            Compute(mesh, *pins_);
        } else if (legacyPinIndices_) {
            ComputeImpl(mesh,
                        autoPlacePair(mesh, legacyPinIndices_->first, legacyPinIndices_->second));
        } else {
            Compute(mesh);
        }
    }

    /**
     * @brief Compute the parameterized mesh using automatic pin selection
     *
     * Selects the first boundary vertex and its boundary-edge neighbor as
     * pinned vertices, placing pin0 at the UV origin and pin1 at distance
     * `|p1 - p0|` along the dominant world-axis of `(p1 - p0)`.
     *
     * @throws MeshException If pin selection fails (no boundary vertices).
     * @throws SolverException If matrix cannot be decomposed or if solver fails
     * to find a solution.
     */
    static void Compute(typename Mesh::Pointer& mesh) { ComputeImpl(mesh, autoSelectPins(mesh)); }

    /**
     * @brief Compute the parameterized mesh with caller-specified pin UVs
     *
     * @param mesh Triangle mesh whose vertex positions will be overwritten with
     * computed 2D UV coordinates (z component set to 0).
     * @param pins Sequence of (vertex index, target UV) pairs. Must contain at
     * least two unique, in-range vertex indices. Each pinned vertex's final UV
     * equals the supplied target verbatim.
     * @throws std::invalid_argument If `pins` has fewer than two entries, a
     * duplicate index, or an out-of-range index.
     * @throws SolverException If matrix cannot be decomposed or if solver fails
     * to find a solution.
     */
    static void Compute(typename Mesh::Pointer& mesh, const PinMap& pins)
    {
        validatePins(mesh, pins);
        ComputeImpl(mesh, pins);
    }

    /**
     * @brief Deprecated: compute with an explicit pin pair by index.
     *
     * Builds a 2-entry PinMap using the LSCM axis-snap convention (pin0 at the
     * UV origin, pin1 on the dominant world-axis of `(p1 - p0)`) and
     * dispatches to the PinMap path.
     *
     * @deprecated Prefer `Compute(mesh, PinMap)`. This overload will be
     * removed in version 3.0.
     */
    [[deprecated("Use Compute(mesh, PinMap); will be removed in 3.0")]] static void Compute(
        typename Mesh::Pointer& mesh, std::size_t pin0Idx, std::size_t pin1Idx)
    {
        ComputeImpl(mesh, autoPlacePair(mesh, pin0Idx, pin1Idx));
    }

private:
    /** Optional explicit pin set configured via `setPins()`. */
    std::optional<PinMap> pins_;
    /** Deprecated: legacy two-pin index pair set via `setPinnedVertices`. */
    std::optional<std::pair<std::size_t, std::size_t>> legacyPinIndices_;

    /** Walk the boundary to pick the default pin pair. */
    static auto selectBoundaryPair(const typename Mesh::Pointer& mesh)
        -> std::pair<typename Mesh::VertPtr, typename Mesh::VertPtr>
    {
        auto boundary = mesh->vertices_boundary();
        if (boundary.empty()) {
            throw MeshException("AngleBasedLSCM: mesh has no boundary vertices");
        }
        auto p0 = boundary.front();
        auto e = p0->edge;
        do {
            if (e->pair->is_boundary()) {
                break;
            }
            e = e->pair->next;
        } while (e != p0->edge);
        if (e == p0->edge and not e->pair->is_boundary()) {
            throw MeshException("Pinned vertex not on boundary");
        }
        auto p1 = e->next->vertex;
        return {p0, p1};
    }

    /**
     * @brief Build a 2-entry PinMap from explicit indices using the LSCM
     * axis-snap convention.
     *
     * pin0 at `{0, 0}`, pin1 at signed distance `|p1 - p0|` on whichever world
     * axis the `(p1 - p0)` vector has the largest magnitude.
     */
    static auto autoPlacePair(const typename Mesh::Pointer& mesh, std::size_t p0Idx,
                              std::size_t p1Idx) -> PinMap
    {
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
        return PinMap{{p0Idx, uv0}, {p1Idx, uv1}};
    }

    /**
     * @brief Auto-select two boundary pins and place them via the LSCM
     * axis-snap convention.
     */
    static auto autoSelectPins(const typename Mesh::Pointer& mesh) -> PinMap
    {
        auto [p0, p1] = selectBoundaryPair(mesh);
        return autoPlacePair(mesh, p0->idx, p1->idx);
    }

    /** Validate that the user's PinMap is well-formed for this mesh. */
    static void validatePins(const typename Mesh::Pointer& mesh, const PinMap& pins)
    {
        if (pins.size() < 2) {
            throw std::invalid_argument("AngleBasedLSCM: PinMap requires at least 2 pins");
        }
        auto numVerts = mesh->num_vertices();
        std::unordered_set<std::size_t> seen;
        seen.reserve(pins.size());
        for (const auto& [vIdx, uv] : pins) {
            if (vIdx >= numVerts) {
                throw std::invalid_argument("AngleBasedLSCM: PinMap vertex index out of range");
            }
            if (!seen.insert(vIdx).second) {
                throw std::invalid_argument("AngleBasedLSCM: PinMap has duplicate vertex index");
            }
        }
    }

    /**
     * @brief Core solver: build the LSCM system from the PinMap, solve for free
     * vertices, and write UVs back to the mesh.
     */
    static void ComputeImpl(typename Mesh::Pointer& mesh, const PinMap& pins)
    {
        using SparseMatrix = Eigen::SparseMatrix<T>;
        using DenseMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

        // Pin placement + LSCM system assembly (shared with HierarchicalLSCM)
        auto parts = detail::lscm::buildSystem<T, Mesh>(mesh, pins);

        // Solve for x
        auto x = detail::SolveLeastSquares<SparseMatrix, DenseMatrix, Solver>(parts.A, parts.b);

        // Assign solution to free vertices. Pin vertices are already at their
        // target UVs (written by buildSystem).
        std::unordered_set<std::size_t> pinIdx;
        pinIdx.reserve(pins.size());
        for (const auto& [vIdx, uv] : pins) {
            pinIdx.insert(vIdx);
        }
        for (const auto& v : mesh->vertices()) {
            if (pinIdx.count(v->idx)) {
                continue;
            }
            auto newIdx = 2 * parts.freeIdxTable.at(v->idx);
            v->pos[0] = x(newIdx, 0);
            v->pos[1] = x(newIdx + 1, 0);
            v->pos[2] = T(0);
        }
    }
};

}  // namespace OpenABF
