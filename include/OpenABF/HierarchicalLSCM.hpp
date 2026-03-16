#pragma once

#include <cmath>
#include <optional>
#include <type_traits>
#include <utility>

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCore>

#include "OpenABF/AngleBasedLSCM.hpp"
#include "OpenABF/Exceptions.hpp"
#include "OpenABF/HalfEdgeMesh.hpp"
#include "OpenABF/Math.hpp"

namespace OpenABF
{

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
