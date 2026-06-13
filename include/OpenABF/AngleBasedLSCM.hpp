#pragma once

#include <optional>
#include <type_traits>
#include <utility>

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>

#include "OpenABF/Exceptions.hpp"
#include "OpenABF/HalfEdgeMesh.hpp"
#include "OpenABF/Math.hpp"
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
 * and templated on Eigen::SparseMatrix<T>. The default SparseLU is robust but
 * slow for large meshes. For iterative solving, prefer
 * `Eigen::ConjugateGradient<Eigen::SparseMatrix<T>, Eigen::Lower|Eigen::Upper>`
 * over the default `Lower`-only variant: the `Lower|Upper` template argument
 * enables Eigen's full-matrix SpMV code path, which is faster and — when
 * compiled with OpenMP — multi-threaded. Using only `Lower` (the Eigen
 * default) routes through `selfadjointView<Lower>`, which is a different
 * internal code path that is never OpenMP-parallelized regardless of
 * `Eigen::setNbThreads()`.
 */
template <typename T, class MeshType = HalfEdgeMesh<T>,
          class Solver = Eigen::SparseLU<Eigen::SparseMatrix<T>, Eigen::COLAMDOrdering<int>>,
          std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
class AngleBasedLSCM
{
public:
    /** @brief Mesh type alias */
    using Mesh = MeshType;

    /** @brief Set the pinned vertex indices used by compute() */
    void setPinnedVertices(std::size_t pin0Idx, std::size_t pin1Idx)
    {
        pinnedVertices_ = {pin0Idx, pin1Idx};
    }

    /** @copydoc AngleBasedLSCM::Compute() */
    void compute(typename Mesh::Pointer& mesh) const
    {
        if (pinnedVertices_) {
            Compute(mesh, pinnedVertices_->first, pinnedVertices_->second);
        } else {
            Compute(mesh);
        }
    }

    /**
     * @brief Compute the parameterized mesh using automatic pin selection
     *
     * Selects the first boundary vertex and its boundary-edge neighbor as
     * pinned vertices.
     *
     * @throws MeshException If pinned vertex is not on boundary.
     * @throws SolverException If matrix cannot be decomposed or if solver fails
     * to find a solution.
     */
    static void Compute(typename Mesh::Pointer& mesh)
    {
        // Pinned vertex selection: first boundary vertex + boundary-edge neighbor
        auto p0 = mesh->vertices_boundary().front();
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
        ComputeImpl(mesh, p0, p1);
    }

    /**
     * @brief Compute the parameterized mesh with explicit pinned vertex indices
     *
     * @param mesh Triangle mesh whose vertex positions will be overwritten with
     * computed 2D UV coordinates (z component set to 0).
     * @param pin0Idx Index of the first pinned vertex (placed at the UV origin)
     * @param pin1Idx Index of the second pinned vertex (placed on the nearest axis)
     * @throws SolverException If matrix cannot be decomposed or if solver fails
     * to find a solution.
     */
    static void Compute(typename Mesh::Pointer& mesh, std::size_t pin0Idx, std::size_t pin1Idx)
    {
        ComputeImpl(mesh, mesh->vertex(pin0Idx), mesh->vertex(pin1Idx));
    }

private:
    /** Optional explicit pin pair set via setPinnedVertices() */
    std::optional<std::pair<std::size_t, std::size_t>> pinnedVertices_;

    /**
     * @brief Core solver: place p0/p1 on the UV axes then solve for free vertices
     */
    static void ComputeImpl(typename Mesh::Pointer& mesh, const typename Mesh::VertPtr& p0,
                            const typename Mesh::VertPtr& p1)
    {
        using SparseMatrix = Eigen::SparseMatrix<T>;
        using DenseMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

        // Pin placement + LSCM system assembly (shared with HierarchicalLSCM)
        auto parts = detail::lscm::buildSystem<T, Mesh>(mesh, p0, p1);

        // Solve for x
        auto x = detail::SolveLeastSquares<SparseMatrix, DenseMatrix, Solver>(parts.A, parts.b);

        // Assign solution to UV coordinates
        // Pins are already updated by buildSystem, so these are free vertices
        for (const auto& v : mesh->vertices()) {
            if (v == p0 or v == p1) {
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