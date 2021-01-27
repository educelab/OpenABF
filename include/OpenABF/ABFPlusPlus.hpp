#pragma once

#include <cmath>

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
 */
template <
    typename T,
    class MeshType = HalfEdgeMesh<
        T,
        3,
        traits::ABFVertexTraits<T>,
        traits::ABFEdgeTraits<T>,
        traits::ABFFaceTraits<T>>,
    std::enable_if_t<std::is_floating_point<T>::value, bool> = true>
class ABFPlusPlus
{
public:
    /** @brief Mesh type alias */
    using Mesh = MeshType;

    /** @brief Set the mesh to be processed */
    void setMesh(typename Mesh::Pointer m) { mesh_ = m; }

    /** @brief Set the maximum number of iterations */
    void setMaxIterations(std::size_t it) { maxIters_ = it; }

    /**
     * @brief Get the mesh gradient
     *
     * **Note:** Result is only valid after running compute().
     */
    auto gradient() const -> T { return grad_; }

    /**
     * @brief Get the number of iterations of the last computation
     *
     * **Note:** Result is only valid after running compute().
     */
    auto iterations() const -> std::size_t { return iters_; }

    /** @brief Compute the interior angles */
    void compute()
    {
        using namespace detail::ABF;
        // TODO: Deep copy mesh

        // Initialize and bound angle properties
        static constexpr auto MinAngle = PI<T> / T(180);
        static constexpr auto MaxAngle = PI<T> - MinAngle;
        for (auto& e : mesh_->edges()) {
            e->alpha = e->beta = e->phi =
                std::min(std::max(e->alpha, MinAngle), MaxAngle);
            e->alpha_sin = std::sin(e->alpha);
            e->alpha_cos = std::cos(e->alpha);
            e->weight = T(1) / (e->phi * e->phi);
        }

        // Update weights for interior vertices
        for (auto& v : mesh_->vertices_interior()) {
            auto wheel = v->wheel();
            auto angle_sum = std::accumulate(
                wheel.begin(), wheel.end(), T(0),
                [](auto a, auto b) { return a + b->beta; });
            for (auto& e : wheel) {
                e->phi *= 2 * PI<T> / angle_sum;
                e->weight = T(1) / (e->phi * e->phi);
            }
        }

        // while ||∇F(x)|| > ε
        grad_ = Gradient<T>(mesh_);
        auto gradDelta = INF<T>;
        iters_ = 0;
        while (grad_ > 0.001 and gradDelta > 0.001 and iters_ < maxIters_) {
            // Typedefs
            using Triplet = Eigen::Triplet<T>;
            using SparseMatrix = Eigen::SparseMatrix<T>;
            using DenseVector = Eigen::Matrix<T, Eigen::Dynamic, 1>;
            using Ordering = Eigen::COLAMDOrdering<int>;
            using Solver = Eigen::SparseLU<SparseMatrix, Ordering>;

            // Helpful parameters
            auto vCnt = mesh_->num_vertices();
            auto vIntCnt = mesh_->num_vertices_interior();
            auto edgeCnt = mesh_->num_edges();
            auto faceCnt = mesh_->num_faces();

            // b1 = -alpha gradient
            std::vector<Triplet> triplets;
            std::size_t idx{0};
            for (const auto& e : mesh_->edges()) {
                triplets.emplace_back(idx, 0, -AlphaGrad<T>(e));
                ++idx;
            }
            SparseMatrix b1(edgeCnt, 1);
            b1.reserve(triplets.size());
            b1.setFromTriplets(triplets.begin(), triplets.end());

            // b2 = -lambda gradient
            triplets.clear();
            idx = 0;
            // lambda tri
            for (const auto& f : mesh_->faces()) {
                triplets.emplace_back(idx, 0, -TriGrad<T>(f));
                idx++;
            }
            // lambda plan and lambda len
            for (const auto& v : mesh_->vertices_interior()) {
                triplets.emplace_back(idx, 0, -PlanGrad<T>(v));
                triplets.emplace_back(vIntCnt + idx, 0, -LenGrad<T>(v));
                idx++;
            }
            SparseMatrix b2(faceCnt + 2 * vIntCnt, 1);
            b2.reserve(triplets.size());
            b2.setFromTriplets(triplets.begin(), triplets.end());

            // vertex idx -> interior vertex idx permutation
            std::map<std::size_t, std::size_t> vIdx2vIntIdx;
            std::size_t newIdx{0};
            for (const auto& v : mesh_->vertices_interior()) {
                vIdx2vIntIdx[v->idx] = newIdx++;
            }

            // Compute J1 + J2
            triplets.clear();
            idx = 0;
            // Jacobian of the CTri constraints
            for (; idx < faceCnt; idx++) {
                triplets.emplace_back(idx, 3 * idx, 1);
                triplets.emplace_back(idx, 3 * idx + 1, 1);
                triplets.emplace_back(idx, 3 * idx + 2, 1);
            }
            for (const auto& v : mesh_->vertices_interior()) {
                for (const auto& e0 : v->wheel()) {
                    // Jacobian of the CPlan constraint
                    triplets.emplace_back(idx, e0->idx, 1);

                    // Jacobian of the CLen constraint
                    auto e1 = e0->next;
                    auto e2 = e1->next;
                    auto d1 = LenGrad<T>(v, e1);
                    auto d2 = LenGrad<T>(v, e2);
                    triplets.emplace_back(vIntCnt + idx, e1->idx, d1);
                    triplets.emplace_back(vIntCnt + idx, e2->idx, d2);
                }
                ++idx;
            }
            SparseMatrix J(faceCnt + 2 * vIntCnt, 3 * faceCnt);
            J.reserve(triplets.size());
            J.setFromTriplets(triplets.begin(), triplets.end());

            // Lambda = diag(2/w)
            // v.weight == 1/w, so LambdaInv is diag(2*weight)
            // We only need Lambda Inverse, so this is 1 / 2*weight
            triplets.clear();
            idx = 0;
            for (const auto& e : mesh_->edges()) {
                triplets.emplace_back(idx, idx, T(1) / (2 * e->weight));
                ++idx;
            }
            SparseMatrix LambdaInv(edgeCnt, edgeCnt);
            LambdaInv.reserve(edgeCnt);
            LambdaInv.setFromTriplets(triplets.begin(), triplets.end());

            // solve Eq. 16
            auto bstar = J * LambdaInv * b1 - b2;
            auto JLiJt = J * LambdaInv * J.transpose();

            SparseMatrix LambdaStarInv = JLiJt.block(0, 0, faceCnt, faceCnt);
            for (std::size_t k=0; k < LambdaStarInv.outerSize(); ++k) {
                for (typename SparseMatrix::InnerIterator it(LambdaStarInv, k);
                     it; ++it) {
                    it.valueRef() = 1.F / it.value();
                }
            }
            auto Jstar = JLiJt.block(faceCnt,0,2*vIntCnt,faceCnt);
            auto JstarT = JLiJt.block(0,faceCnt,faceCnt, 2*vIntCnt);
            auto Jstar2 = JLiJt.block(faceCnt,faceCnt,2*vIntCnt, 2*vIntCnt);
            auto bstar1 = bstar.block(0, 0, faceCnt, 1);
            auto bstar2 = bstar.block(faceCnt, 0, 2*vIntCnt, 1);

            // (J* Lam*^-1 J*^t - J**) delta_lambda_2 = J* Lam*^-1 b*_1 - b*_2
            SparseMatrix A = Jstar * LambdaStarInv * JstarT - Jstar2;
            SparseMatrix b = Jstar * LambdaStarInv * bstar1 - bstar2;
            A.makeCompressed();
            Solver solver;
            solver.compute(A);
            if (solver.info() != Eigen::ComputationInfo::Success) {
                throw SolverException(solver.lastErrorMessage());
            }
            auto deltaLambda2 = solver.solve(b);
            if (solver.info() != Eigen::ComputationInfo::Success) {
                throw SolverException(solver.lastErrorMessage());
            }

            // Compute Eq. 17 -> delta_lambda_1
            auto deltaLambda1 =
                LambdaStarInv * (bstar1 - JstarT * deltaLambda2);

            // Construct deltaLambda
            DenseVector deltaLambda(
                deltaLambda1.rows() + deltaLambda2.rows(), 1);
            deltaLambda << DenseVector(deltaLambda1), DenseVector(deltaLambda2);

            // Compute Eq. 10 -> delta_alpha
            DenseVector deltaAlpha =
                LambdaInv * (b1 - J.transpose() * deltaLambda);

            // lambda += delta_lambda
            for (auto& f : mesh_->faces()) {
                f->lambda_tri += deltaLambda(f->idx, 0);
            }
            for (auto& v : mesh_->vertices_interior()) {
                auto intIdx = vIdx2vIntIdx.at(v->idx);
                v->lambda_plan += deltaLambda(faceCnt + intIdx, 0);
                v->lambda_len += deltaLambda(faceCnt + vIntCnt + intIdx, 0);
            }

            // alpha += delta_alpha
            // Update sin and cos
            idx = 0;
            for (auto& e : mesh_->edges()) {
                e->alpha += deltaAlpha(idx++, 0);
                e->alpha = std::min(std::max(e->alpha, T(0)), PI<T>);
                e->alpha_sin = std::sin(e->alpha);
                e->alpha_cos = std::cos(e->alpha);
            }

            // Recalculate gradient for next iteration
            auto newGrad = Gradient<T>(mesh_);
            gradDelta = std::abs(newGrad - grad_);
            grad_ = newGrad;
            iters_++;
        }
    }
private:
    /** Mesh */
    typename Mesh::Pointer mesh_;
    /** Gradient */
    T grad_{0};
    /** Number of executed iterations */
    std::size_t iters_{0};
    /** Max iterations */
    std::size_t maxIters_{10};
};

}  // namespace OpenABF