#pragma once

#include <cmath>
#include <map>

#include <Eigen/SparseLU>

#include "OpenABF/HalfEdgeMesh.hpp"
#include "OpenABF/Math.hpp"

namespace OpenABF
{

/**
 * @brief Compute parameterized mesh using Angle-based LSCM
 *
 * Computes a least-squares conformal parameterization of a mesh. Unlike the
 * original LSCM algorithm, this class ignores the vertex positions and instead
 * uses the associated edge angles (MeshType::EdgeTraits::alpha) to compute the
 * initial vertex positions. The parameterization can be improved by processing
 * the mesh with a parameterized angle optimizer, such as ABFPlusPlus, before
 * processing with this class.
 *
 * Implements the angle-based variant of "Least squares conformal maps for
 * automatic texture atlas generation" by Lévy _et al._ (2002)
 * \cite levy2002lscm.
 *
 * @tparam T
 * @tparam MeshType
 */
template <
    typename T,
    class MeshType = HalfEdgeMesh<
        T,
        3,
        traits::DefaultVertexTraits<T>,
        traits::DefaultEdgeTraits<T>,
        traits::DefaultFaceTraits<T>>,
    std::enable_if_t<std::is_floating_point<T>::value, bool> = true>
class AngleBasedLSCM
{
public:
    /** @brief Mesh type alias */
    using Mesh = MeshType;

    /** @brief Set the mesh to be processed */
    void setMesh(const typename Mesh::Pointer& m) { mesh_ = m; }

    /** @brief Compute the parameterized mesh */
    void compute()
    {
        using Triplet = Eigen::Triplet<T>;
        using SparseMatrix = Eigen::SparseMatrix<T>;
        using DenseMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
        using Ordering = Eigen::COLAMDOrdering<int>;
        using Solver = Eigen::SparseLU<SparseMatrix, Ordering>;

        // Pinned vertex selection
        // Get the end points of a boundary edge
        // TODO: Initial pin selection could be better
        auto p0 = mesh_->vertices_boundary()[0];
        auto e = p0->edge;
        do {
            if (not e->pair) {
                break;
            }
            e = e->pair->next;
        } while (e != p0->edge);
        if (e == p0->edge) {
            throw std::invalid_argument("Vertex not on boundary");
        }
        auto p1 = e->next->vertex;

        // Map selected edge to closest XY axis
        // TODO: Project this triangle's plane to UV?
        auto pinVec = p1->pos - p0->pos;
        auto dist = norm(pinVec);
        pinVec /= dist;
        p0->pos = {T(0), T(0), T(0)};
        auto maxElem = std::max_element(pinVec.begin(), pinVec.end());
        auto maxAxis = std::distance(pinVec.begin(), maxElem);
        if (maxAxis == 0) {
            p1->pos = {dist, T(0), T(0)};
        } else {
            p1->pos = {T(0), dist, T(0)};
        }

        // For convenience
        auto numFaces = mesh_->num_faces();
        auto numVerts = mesh_->num_vertices();
        auto numFixed = 2;
        auto numFree = numVerts - numFixed;

        // Permutation for free vertices
        // This helps us find a vert's row in the solution matrix
        std::map<std::size_t, std::size_t> freeIdxTable;
        for (const auto& v : mesh_->vertices()) {
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
        for (const auto& f : mesh_->faces()) {
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

            // Rotate the edges of the triangle so the final one has max sin
            // Blender's implementation notes that this is stable ordering
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

            auto row = f->idx;
            // If pin0 or pin1, put in fixedB matrix, else put in A
            // TODO: This can be improved
            if (e0->vertex == p0) {
                tripletsB.emplace_back(row, 0, cosine - T(1));
                tripletsB.emplace_back(row, 1, -sine);
                tripletsB.emplace_back(row + 1, 0, sine);
                tripletsB.emplace_back(row + 1, 1, cosine - T(1));
            } else if (e0->vertex == p1) {
                tripletsB.emplace_back(row, 2, cosine - T(1));
                tripletsB.emplace_back(row, 3, -sine);
                tripletsB.emplace_back(row + 1, 2, sine);
                tripletsB.emplace_back(row + 1, 3, cosine - T(1));
            } else {
                auto freeIdx = freeIdxTable.at(e0->vertex->idx);
                tripletsA.emplace_back(row, 2 * freeIdx, cosine - T(1));
                tripletsA.emplace_back(row, 2 * freeIdx + 1, -sine);
                tripletsA.emplace_back(row + 1, 2 * freeIdx, sine);
                tripletsA.emplace_back(row + 1, 2 * freeIdx + 1, cosine - T(1));
            }

            if (e1->vertex == p0) {
                tripletsB.emplace_back(row, 0, -cosine);
                tripletsB.emplace_back(row, 1, sine);
                tripletsB.emplace_back(row + 1, 0, -sine);
                tripletsB.emplace_back(row + 1, 1, -cosine);
            } else if (e1->vertex == p1) {
                tripletsB.emplace_back(row, 2, -cosine);
                tripletsB.emplace_back(row, 3, sine);
                tripletsB.emplace_back(row + 1, 2, -sine);
                tripletsB.emplace_back(row + 1, 3, -cosine);
            } else {
                auto freeIdx = freeIdxTable.at(e1->vertex->idx);
                tripletsA.emplace_back(row, 2 * freeIdx, -cosine);
                tripletsA.emplace_back(row, 2 * freeIdx + 1, sine);
                tripletsA.emplace_back(row + 1, 2 * freeIdx, -sine);
                tripletsA.emplace_back(row + 1, 2 * freeIdx + 1, -cosine);
            }

            if (e2->vertex == p0) {
                tripletsB.emplace_back(row, 0, T(1));
                tripletsB.emplace_back(row + 1, 1, T(1));
            } else if (e2->vertex == p1) {
                tripletsB.emplace_back(row, 2, T(1));
                tripletsB.emplace_back(row + 1, 3, T(1));
            } else {
                auto freeIdx = freeIdxTable.at(e2->vertex->idx);
                tripletsA.emplace_back(row, 2 * freeIdx, T(1));
                tripletsA.emplace_back(row + 1, 2 * freeIdx + 1, T(1));
            }
        }
        SparseMatrix A(2 * numFaces, 2 * numFree);
        A.reserve(tripletsA.size());
        A.setFromTriplets(tripletsA.begin(), tripletsA.end());

        SparseMatrix bFree(2 * numFaces, 2 * numFixed);
        bFree.reserve(tripletsB.size());
        bFree.setFromTriplets(tripletsB.begin(), tripletsB.end());

        // Calculate rhs from free and fixed matrices
        SparseMatrix b = bFree * bFixed * -1;

        // Setup AtA and solver
        SparseMatrix AtA = A.transpose() * A;
        AtA.makeCompressed();
        Solver solver;
        solver.compute(AtA);
        if (solver.info() != Eigen::ComputationInfo::Success) {
            throw SolverException(solver.lastErrorMessage());
        }

        // Setup Atb
        SparseMatrix Atb = A.transpose() * b;

        // Solve AtAx = AtAb
        DenseMatrix x = solver.solve(Atb);

        // Assign solution to UV coordinates
        // Pins are already updated, so these are free vertices
        for (const auto& v : mesh_->vertices()) {
            if (v == p0 or v == p1) {
                continue;
            }
            auto newIdx = 2 * freeIdxTable.at(v->idx);
            v->pos[0] = x(newIdx, 0);
            v->pos[1] = x(newIdx + 1, 0);
            v->pos[2] = T(0);
        }
    }

private:
    /** Mesh */
    typename Mesh::Pointer mesh_;
};

}  // namespace OpenABF