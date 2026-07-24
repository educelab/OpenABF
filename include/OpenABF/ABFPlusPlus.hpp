#pragma once

#include <cassert>
#include <cmath>
#include <limits>
#include <vector>

#include <Eigen/SparseLU>

#include "OpenABF/ABF.hpp"
#include "OpenABF/Exceptions.hpp"
#include "OpenABF/HalfEdgeMesh.hpp"
#include "OpenABF/Math.hpp"

namespace OpenABF
{

/**
 * @brief Compute parameterized interior angles using ABF++
 *
 * Iteratively computes a new set of interior angles which minimize the total
 * angular error of the parameterized mesh. This follows the ABF++ formulation,
 * which solves a 5x smaller system of equations than standard ABF at the
 * expense of more iterations.
 *
 * This class **does not** compute a parameterized mesh. Rather, it calculates
 * the optimal interior angles for such a mesh. To convert this information
 * into a full parameterization, pass the processed HalfEdgeMesh to
 * AngleBasedLSCM.
 *
 * Implements "ABF++: Fast and Robust Angle Based Flattening" by Sheffer
 * _et al._ (2005) \cite sheffer2005abf++.
 *
 * @tparam T Floating-point type
 * @tparam MeshType HalfEdgeMesh type which implements the ABF traits
 * @tparam Solver A solver implementing the
 * [Eigen Sparse solver
 * concept](https://eigen.tuxfamily.org/dox-devel/group__TopicSparseSystems.html)
 * and templated on Eigen::SparseMatrix<T>
 */
template <typename T, class MeshType = detail::ABF::Mesh<T>,
          class Solver = Eigen::SparseLU<Eigen::SparseMatrix<T>, Eigen::COLAMDOrdering<int>>,
          std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
class ABFPlusPlus
{
public:
    /** @brief Mesh type alias */
    using Mesh = MeshType;

    /** @brief Set the maximum number of iterations */
    void set_max_iterations(std::size_t it) { max_iters_ = it; }

    /** @brief Set the gradient convergence threshold */
    void set_gradient_threshold(T t) { grad_threshold_ = t; }

    /** @deprecated Use `set_max_iterations`; will be removed in 3.0. */
    [[deprecated("Use set_max_iterations; will be removed in 3.0")]] void setMaxIterations(
        std::size_t it)
    {
        set_max_iterations(it);
    }

    /** @deprecated Use `set_gradient_threshold`; will be removed in 3.0. */
    [[deprecated("Use set_gradient_threshold; will be removed in 3.0")]] void setGradientThreshold(
        T t)
    {
        set_gradient_threshold(t);
    }

    /**
     * @brief Get the mesh gradient
     *
     * **Note:** Result is only valid after running compute().
     */
    [[nodiscard]] auto gradient() const -> T { return grad_; }

    /**
     * @brief Get the number of iterations of the last computation
     *
     * **Note:** Result is only valid after running compute().
     */
    [[nodiscard]] auto iterations() const -> std::size_t { return iters_; }

    /** @copydoc ABFPlusPlus::Compute */
    void compute(typename Mesh::Pointer& mesh)
    {
        Compute(mesh, iters_, grad_, max_iters_, grad_threshold_);
    }

    /**
     * @brief Compute parameterized interior angles
     *
     * @throws SolverException If matrix cannot be decomposed or if solver fails
     * to find a solution.
     * @throws MeshException If mesh gradient cannot be calculated.
     */
    static void Compute(typename Mesh::Pointer& mesh, std::size_t& iters, T& gradient,
                        const std::size_t maxIters = 10, T gradThreshold = T(0.001))
    {
        using namespace detail::ABF;

        // Initialize angles and weights
        InitializeAnglesAndWeights<T>(mesh);

        // while ||∇F(x)|| > ε
        gradient = Gradient<T>(mesh);
        if (std::isnan(gradient) or std::isinf(gradient)) {
            throw MeshException("Mesh gradient cannot be computed");
        }
        auto gradDelta = INF<T>;
        iters = 0;

        // vertex idx -> interior vertex idx lookup (pre-built once, O(1) access)
        std::vector<std::size_t> vIdx2vIntIdx(mesh->num_vertices(),
                                              std::numeric_limits<std::size_t>::max());
        {
            std::size_t newIdx{0};
            for (const auto& v : mesh->vertices_interior()) {
                vIdx2vIntIdx[v->idx] = newIdx++;
            }
        }

        // Typedefs
        using Triplet = Eigen::Triplet<T>;
        using SparseMatrix = Eigen::SparseMatrix<T>;
        using DenseVector = Eigen::Matrix<T, Eigen::Dynamic, 1>;

        // Helpful parameters
        auto vIntCnt = mesh->num_vertices_interior();
        auto edgeCnt = mesh->num_edges();
        auto faceCnt = mesh->num_faces();

        std::vector<Triplet> triplets;
        SparseMatrix b1(edgeCnt, 1);
        SparseMatrix b2(faceCnt + 2 * vIntCnt, 1);
        SparseMatrix J(faceCnt + 2 * vIntCnt, 3 * faceCnt);
        SparseMatrix LambdaInv(edgeCnt, edgeCnt);
        SparseMatrix LambdaStarInv(faceCnt, faceCnt);
        SparseMatrix A(2 * vIntCnt, 2 * vIntCnt);
        SparseMatrix b(2 * vIntCnt, 1);
        DenseVector deltaLambda(faceCnt + 2 * vIntCnt, 1);
        DenseVector deltaAlpha(edgeCnt, 1);

        while (gradient > gradThreshold and gradDelta > gradThreshold and iters < maxIters) {
            if (std::isnan(gradient) or std::isinf(gradient)) {
                throw MeshException("Mesh gradient cannot be computed");
            }

            // b1 = -alpha gradient
            triplets.clear();
            std::size_t idx{0};
            for (const auto& e : mesh->edges()) {
                triplets.emplace_back(idx, 0, -AlphaGrad<T>(e));
                ++idx;
            }
            b1.setFromTriplets(triplets.begin(), triplets.end());

            // b2 = -lambda gradient
            triplets.clear();
            idx = 0;
            // lambda tri
            for (const auto& f : mesh->faces()) {
                triplets.emplace_back(idx, 0, -TriGrad<T>(f));
                idx++;
            }
            // lambda plan and lambda len
            for (const auto& v : mesh->vertices_interior()) {
                triplets.emplace_back(idx, 0, -PlanGrad<T>(v));
                triplets.emplace_back(vIntCnt + idx, 0, -LenGrad<T>(v));
                idx++;
            }
            b2.setFromTriplets(triplets.begin(), triplets.end());

            // Compute J1 + J2
            triplets.clear();
            idx = 0;
            // Jacobian of the CTri constraints
            for (; idx < faceCnt; idx++) {
                triplets.emplace_back(idx, 3 * idx, 1);
                triplets.emplace_back(idx, 3 * idx + 1, 1);
                triplets.emplace_back(idx, 3 * idx + 2, 1);
            }
            for (const auto& v : mesh->vertices_interior()) {
                for (const auto& e0 : v->wheel()) {
                    // Jacobian of the CPlan constraint
                    triplets.emplace_back(idx, e0->idxI.value(), 1);

                    // Jacobian of the CLen constraint
                    auto e1 = e0->next;
                    auto e2 = e1->next;
                    auto d1 = LenGrad<T>(v, e1);
                    auto d2 = LenGrad<T>(v, e2);
                    triplets.emplace_back(vIntCnt + idx, e1->idxI.value(), d1);
                    triplets.emplace_back(vIntCnt + idx, e2->idxI.value(), d2);
                }
                ++idx;
            }
            J.setFromTriplets(triplets.begin(), triplets.end());

            // Lambda = diag(2/w)
            // v.weight == 1/w, so LambdaInv is diag(2*weight)
            // We only need Lambda Inverse, so this is 1 / 2*weight
            triplets.clear();
            idx = 0;
            for (const auto& e : mesh->edges()) {
                triplets.emplace_back(idx, idx, T(1) / (2 * e->weight));
                ++idx;
            }
            LambdaInv.setFromTriplets(triplets.begin(), triplets.end());

            // solve Eq. 16
            auto bstar = J * LambdaInv * b1 - b2;
            auto JLiJt = J * LambdaInv * J.transpose();

            LambdaStarInv = JLiJt.block(0, 0, faceCnt, faceCnt);
            for (int k = 0; k < LambdaStarInv.outerSize(); ++k) {
                for (typename SparseMatrix::InnerIterator it(LambdaStarInv, k); it; ++it) {
                    it.valueRef() = T(1) / it.value();
                }
            }
            auto Jstar = JLiJt.block(faceCnt, 0, 2 * vIntCnt, faceCnt);
            auto JstarT = JLiJt.block(0, faceCnt, faceCnt, 2 * vIntCnt);
            auto Jstar2 = JLiJt.block(faceCnt, faceCnt, 2 * vIntCnt, 2 * vIntCnt);
            auto bstar1 = bstar.block(0, 0, faceCnt, 1);
            auto bstar2 = bstar.block(faceCnt, 0, 2 * vIntCnt, 1);

            // (J* Lam*^-1 J*^t - J**) delta_lambda_2 = J* Lam*^-1 b*_1 - b*_2
            A = Jstar * LambdaStarInv * JstarT - Jstar2;
            b = Jstar * LambdaStarInv * bstar1 - bstar2;
            A.makeCompressed();
            Solver solver;
            solver.compute(A);
            if (solver.info() != Eigen::ComputationInfo::Success) {
                throw SolverException("ABF++: Failed to solve A");
            }
            auto deltaLambda2 = solver.solve(b);
            if (solver.info() != Eigen::ComputationInfo::Success) {
                throw SolverException("ABF++: Failed to solve b");
            }

            // Compute Eq. 17 -> delta_lambda_1
            auto deltaLambda1 = LambdaStarInv * (bstar1 - JstarT * deltaLambda2);

            // Construct deltaLambda
            deltaLambda.topRows(faceCnt) = deltaLambda1;
            deltaLambda.bottomRows(2 * vIntCnt) = deltaLambda2;

            // Compute Eq. 10 -> delta_alpha
            deltaAlpha = LambdaInv * (b1 - J.transpose() * deltaLambda);

            // lambda += delta_lambda
            for (auto& f : mesh->faces()) {
                f->lambda_tri += deltaLambda(f->idx, 0);
            }
            for (auto& v : mesh->vertices_interior()) {
                auto intIdx = vIdx2vIntIdx[v->idx];
                assert(intIdx != std::numeric_limits<std::size_t>::max());
                v->lambda_plan += deltaLambda(faceCnt + intIdx, 0);
                v->lambda_len += deltaLambda(faceCnt + vIntCnt + intIdx, 0);
            }

            // alpha += delta_alpha
            // Update sin and cos
            idx = 0;
            for (auto& e : mesh->edges()) {
                e->alpha += deltaAlpha(idx++, 0);
                e->alpha = std::min(std::max(e->alpha, T(0)), PI<T>);
                e->alpha_sin = std::sin(e->alpha);
                e->alpha_cos = std::cos(e->alpha);
            }

            // Recalculate gradient for next iteration
            auto newGrad = Gradient<T>(mesh);
            gradDelta = std::abs(newGrad - gradient);
            gradient = newGrad;
            iters++;
        }
    }

    /** @brief Compute parameterized interior angles */
    static void Compute(typename Mesh::Pointer& mesh)
    {
        std::size_t iters{0};
        T grad{0};
        Compute(mesh, iters, grad);
    }

private:
    /** Gradient */
    T grad_{0};
    /** Number of executed iterations */
    std::size_t iters_{0};
    /** Max iterations */
    std::size_t max_iters_{10};
    /** Gradient convergence threshold */
    T grad_threshold_{0.001};
};

}  // namespace OpenABF