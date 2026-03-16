#pragma once

#include <cmath>
#include <map>
#include <optional>
#include <type_traits>
#include <utility>

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>

#include "OpenABF/Exceptions.hpp"
#include "OpenABF/HalfEdgeMesh.hpp"
#include "OpenABF/Math.hpp"

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
 * and templated on Eigen::SparseMatrix<T>
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
        using Triplet = Eigen::Triplet<T>;
        using SparseMatrix = Eigen::SparseMatrix<T>;
        using DenseMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

        // Map selected edge to closest XY axis
        // Use sign to select direction
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

        // For convenience
        auto numFaces = mesh->num_faces();
        auto numVerts = mesh->num_vertices();
        auto numFixed = 2;
        auto numFree = numVerts - numFixed;

        // Permutation for free vertices
        // This helps us find a vert's row in the solution matrix
        std::map<std::size_t, std::size_t> freeIdxTable;
        for (const auto& v : mesh->vertices()) {
            if (v == p0 or v == p1) {
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

        // Setup variables matrix
        // Are only solving for free vertices, so push pins in special matrix
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

            // Find the max sin idx
            std::vector<T> sins{sin0, sin1, sin2};
            auto sinMaxElem = std::max_element(sins.begin(), sins.end());
            auto sinMaxIdx = std::distance(sins.begin(), sinMaxElem);

            // Rotate the edge order of the face so last angle is largest
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

            // Assemble per-vertex contributions for this face (Lévy et al. 2002, Eq. 10)
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

        // Calculate rhs from free and fixed matrices
        SparseMatrix b = bFree * bFixed * -1;

        // Solve for x
        auto x = detail::SolveLeastSquares<SparseMatrix, DenseMatrix, Solver>(A, b);

        // Assign solution to UV coordinates
        // Pins are already updated, so these are free vertices
        for (const auto& v : mesh->vertices()) {
            if (v == p0 or v == p1) {
                continue;
            }
            auto newIdx = 2 * freeIdxTable.at(v->idx);
            v->pos[0] = x(newIdx, 0);
            v->pos[1] = x(newIdx + 1, 0);
            v->pos[2] = T(0);
        }
    }
};

}  // namespace OpenABF