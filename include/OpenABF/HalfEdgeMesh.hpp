#pragma once

#include <algorithm>
#include <array>
#include <iterator>
#include <map>
#include <memory>
#include <vector>

#include "OpenABF/Vector.hpp"

namespace OpenABF
{

namespace traits
{
/** @brief Default HalfEdgeMesh vertex traits */
template <typename T>
struct DefaultVertexTraits {
};

/** @brief Default HalfEdgeMesh edge traits */
template <typename T>
struct DefaultEdgeTraits {
    /** Edge internal angle */
    T alpha{0};
};

/** @brief Default HalfEdgeMesh face traits */
template <typename T>
struct DefaultFaceTraits {
};
}  // namespace traits

/**
 * @brief Compute the internal angles of a face
 *
 * Updates the current angle (DefaultEdgeTraits::alpha) with the internal angles
 * derived from the face's vertex positions. Useful if you want to reset a face
 * after being processed by ABF or ABFPlusPlus.
 *
 * @tparam FacePtr A Face-type pointer implementing DefaultEdgeTraits
 */
template <class FacePtr>
void ComputeFaceAngles(FacePtr& face)
{
    for (auto& e : *face) {
        auto ab = e->next->vertex->pos - e->vertex->pos;
        auto ac = e->next->next->vertex->pos - e->vertex->pos;
        e->alpha = interior_angle(ab, ac);
    }
}

/**
 * @brief Compute the internal angles for all faces in a mesh
 *
 * Runs ComputeFaceAngles on all faces in the mesh. Useful if you want to reset
 * a mesh after running through ABF or ABFPlusPlus.
 *
 * @tparam MeshPtr A Mesh-type pointer with faces implementing DefaultEdgeTraits
 */
template <class MeshPtr>
void ComputeMeshAngles(MeshPtr& mesh)
{
    for (auto& f : mesh->faces()) {
        ComputeFaceAngles(f);
    }
}

/**
 * @brief Half-edge mesh class
 *
 * A half-edge mesh represents each edge as two oppositely oriented half-edges.
 * There is one half-edge for each face containing the original edge. For
 * example, if two faces share the edge **AB**, this will result in two
 * half-edges, **AB** and **BA**. If an edge **BC** lies on the mesh border
 * (i.e. it is only included in a single face), there will one be a single
 * half-edge created. This data structure makes it possible to easily traverse
 * the edges and faces adjacent to a vertex (the "wheel"), as well as to
 * traverse the edges of each face.
 *
 * For more information, see Chapter 12.1 in "Fundamentals of Computer
 * Graphics, Fourth edition", Marschner and Shirley (2015)
 * \cite marschner2015fundamentals.
 *
 * @tparam T Floating-point type for coordinates
 * @tparam Dim Dimensionality of vertex coordinates
 * @tparam VertexTraits Additional traits for vertices
 * @tparam EdgeTraits Additional traits for edges
 * @tparam FaceTraits Additional traits for face
 */
template <
    typename T,
    std::size_t Dim = 3,
    typename VertexTraits = traits::DefaultVertexTraits<T>,
    typename EdgeTraits = traits::DefaultEdgeTraits<T>,
    typename FaceTraits = traits::DefaultFaceTraits<T>>
class HalfEdgeMesh
{
public:
    /** Pointer type */
    using Pointer = std::shared_ptr<HalfEdgeMesh>;

    struct Vertex;
    struct Edge;
    struct Face;

    /** @brief Vertex pointer type */
    using VertPtr = std::shared_ptr<Vertex>;
    /** @brief Edge pointer type */
    using EdgePtr = std::shared_ptr<Edge>;
    /** @brief Edge pointer type */
    using FacePtr = std::shared_ptr<Face>;

private:
    /**
     * @brief Iterator for the edges of a face
     *
     * @tparam Const If true, is a const iterator
     */
    template <bool Const = false>
    class FaceIterator
    {
    public:
        /** Difference type */
        using difference_type = std::size_t;
        /** Value type */
        using value_type = EdgePtr;
        /** Pointer type */
        using pointer =
            std::conditional_t<Const, value_type const*, value_type*>;
        /** Reference type */
        using reference =
            std::conditional_t<Const, value_type const&, value_type&>;
        /** Iterator category */
        using iterator_category = std::input_iterator_tag;

        /** Default constructor == End iterator */
        FaceIterator() = default;
        /** Construct from head of triangle and current edge */
        explicit FaceIterator(const EdgePtr& head, const EdgePtr& current)
            : head_{head}, current_{current}
        {
        }

        /** Dereference operator */
        template <bool Const_ = Const>
        std::enable_if_t<Const_, reference> operator*() const
        {
            return current_;
        }

        /** Dereference operator */
        template <bool Const_ = Const>
        std::enable_if_t<not Const_, reference> operator*()
        {
            return current_;
        }

        /** Equality operator */
        auto operator==(const FaceIterator& other) const -> bool
        {
            return current_ == other.current_;
        }
        /** Inequality operator */
        auto operator!=(const FaceIterator& other) const -> bool
        {
            return !(*this == other);
        }
        /** Increment operator */
        auto operator++() -> FaceIterator&
        {
            // Already at end
            if (current_ == nullptr) {
                return *this;
            }

            // Get the next edge
            current_ = current_->next;
            // If back at head, done iterating
            if (current_ == head_) {
                current_ = nullptr;
            }
            return *this;
        }
    private:
        /** Pointer to beginning of face */
        EdgePtr head_;
        /** Current edge pointer */
        EdgePtr current_;
    };

public:
    /** @brief %Vertex type */
    struct Vertex : public VertexTraits {
        /** @brief Default constructor */
        Vertex() = default;

        /** @brief Construct from position values */
        template <typename... Args>
        explicit Vertex(Args... args) : pos{args...}
        {
        }

        /** @brief Construct a new Vertex pointer */
        template <typename... Args>
        static auto New(Args&&... args) -> VertPtr
        {
            return std::make_shared<Vertex>(std::forward<Args>(args)...);
        }

        /** @brief Get the edges of a vertex's wheel */
        auto wheel() const -> std::vector<EdgePtr>
        {
            std::vector<EdgePtr> ret;
            auto e = edge;
            do {
                if (not e->pair) {
                    throw std::runtime_error(
                        "Cannot enumerate wheel of boundary vertex.");
                }
                ret.push_back(e);
                e = e->pair->next;
            } while (e != edge);
            return ret;
        }

        /** @brief Returns if vertex is on mesh boundary */
        auto is_boundary() const -> bool
        {
            auto e = edge;
            do {
                if (not e->pair) {
                    return true;
                }
                e = e->pair->next;
            } while (e != edge);
            return false;
        }

        /** @brief Returns if vertex is interior to mesh */
        auto is_interior() const -> bool { return not is_boundary(); }

        /** Insertion index */
        std::size_t idx{0};
        /** Vertex position */
        Vec<T, Dim> pos;
        /**
         * Pointer to _an_ Edge with this vertex as its head. Note that there
         * may be many such vertices.
         */
        EdgePtr edge;
    };

    /** @brief %Edge type */
    struct Edge : public EdgeTraits {
        /** @brief Construct a new Edge pointer */
        template <typename... Args>
        static EdgePtr New(Args&&... args)
        {
            return std::make_shared<Edge>(std::forward<Args>(args)...);
        }

        /** @brief Returns if edge is on mesh boundary */
        auto is_boundary() const -> bool { return pair == nullptr; }

        /**
         * @brief This edge's adjacent half-edge
         *
         * If nullptr, this edge is on the boundary of the mesh (no adjacent
         * face).
         */
        EdgePtr pair;
        /** @brief The next edge in this edge's face */
        EdgePtr next;
        /** @brief The edge's vertex */
        VertPtr vertex;
        /** @brief The face containing this edge */
        FacePtr face;
        /** @brief Insertion index */
        std::size_t idx;
    };

    /** @brief %Face type */
    struct Face : public FaceTraits {
        /** @brief Construct a new Face pointer */
        template <typename... Args>
        static auto New(Args&&... args) -> FacePtr
        {
            return std::make_shared<Face>(std::forward<Args>(args)...);
        }

        /** Face edge iterator type */
        using iterator = FaceIterator<false>;
        /** Face edge const iterator type */
        using const_iterator = FaceIterator<true>;
        /** @brief Returns an iterator over the edges of the face */
        iterator begin() { return iterator{head, head}; }
        /** @brief Returns the end iterator */
        iterator end() { return iterator(); }
        /** @brief Returns an const iterator over the edges of the face */
        const_iterator cbegin() const { return const_iterator{head, head}; }
        /** @brief Returns the const end iterator */
        const_iterator cend() const { return const_iterator(); }

        /** @brief First edge in the face */
        EdgePtr head;
        /** @brief The next face in the mesh */
        FacePtr next;
        /** @brief Insertion index */
        std::size_t idx;
    };

private:
    /** List of vertices */
    std::vector<VertPtr> verts_;
    /** List of faces */
    std::vector<FacePtr> faces_;
    /** List of edges, indexed by the vertex's insertion index */
    std::multimap<std::size_t, EdgePtr> edges_;

public:
    /** @brief Default constructor */
    HalfEdgeMesh() = default;

    /** @brief Destructor deallocating all element pointers */
    ~HalfEdgeMesh()
    {
        // Remove smart pointers from all items
        for (auto& v : verts_) {
            v->edge = nullptr;
        }
        for (auto& e : edges_) {
            e.second->pair = nullptr;
            e.second->next = nullptr;
            e.second->vertex = nullptr;
            e.second->face = nullptr;
        }
        for (auto& f : faces_) {
            f->head = nullptr;
            f->next = nullptr;
        }
        verts_.clear();
        edges_.clear();
        faces_.clear();
    }

    /** @brief Construct a new HalfEdgeMesh pointer */
    template <typename... Args>
    static Pointer New(Args... args)
    {
        return std::make_shared<HalfEdgeMesh>(std::forward<Args>(args)...);
    }

    /**
     * @brief Insert a new vertex
     *
     * Accepts all arguments supported by the Vertex constructor.
     */
    template <typename... Args>
    auto insert_vertex(Args... args) -> std::size_t
    {
        auto vert = Vertex::New(std::forward<Args>(args)...);
        vert->idx = verts_.size();
        verts_.push_back(vert);
        return vert->idx;
    }

    /**
     * @brief Insert a face from an ordered list of Vertex indices
     *
     * Accepts an iterable supporting range-based for loops.
     */
    template <class Vector>
    auto insert_face(const Vector& vector) -> std::size_t
    {
        // Make a new face structure
        auto face = Face::New();

        // Iterate over the vertex indices
        std::size_t prevIdx{0};
        EdgePtr prevEdge;
        for (const auto& idx : vector) {
            // Make a new edge
            auto newEdge = Edge::New();
            newEdge->face = face;

            // Set the head edge for this face
            if (not face->head) {
                face->head = newEdge;
            }

            // Get the vertex by index
            auto vert = verts_.at(idx);
            newEdge->vertex = vert;
            if (not vert->edge) {
                vert->edge = newEdge;
            }

            // If there's a previous edge
            if (prevEdge) {
                // Update the previous edge's successor
                prevEdge->next = newEdge;

                // Try to find a pair for prev edge using this edge's index
                auto pair = find_edge_(idx, prevIdx);
                if (pair) {
                    if (pair->pair) {
                        throw std::invalid_argument(
                            "Resolved edge already paired.");
                    }
                    prevEdge->pair = pair;
                    pair->pair = prevEdge;
                }
            }

            // Store the edge
            newEdge->idx = edges_.size();
            edges_.emplace(idx, newEdge);

            // Update for the next iteration
            prevIdx = idx;
            prevEdge = newEdge;
        }

        // Link back to the beginning
        prevEdge->next = face->head;

        // Try to find a pair for final edge using this edge's index
        auto pair = find_edge_(face->head->vertex->idx, prevIdx);
        if (pair) {
            if (pair->pair) {
                throw std::invalid_argument("Resolved edge already paired.");
            }
            prevEdge->pair = pair;
            pair->pair = prevEdge;
        }

        // Compute angles for edges in face
        ComputeFaceAngles(face);

        // Give this face an idx and link the previous face with this one
        face->idx = faces_.size();
        if (not faces_.empty()) {
            faces_.back()->next = face;
        }
        faces_.emplace_back(face);
        return face->idx;
    }

    /** @brief Insert a new face from an ordered list of Vertex indices */
    template <typename... Args>
    auto insert_face(Args... args) -> std::size_t
    {
        static_assert(sizeof...(args) >= 3, "Faces require >= 3 indices");
        using Tuple = std::tuple<Args...>;
        using ElemT = typename std::tuple_element<0, Tuple>::type;
        return insert_face(std::initializer_list<ElemT>{args...});
    }

    /** @brief Get the list of vertices in insertion order */
    auto vertices() const -> std::vector<VertPtr> { return verts_; }

    /** @brief Get the list of edges in insertion order */
    auto edges() const -> std::vector<EdgePtr>
    {
        std::vector<EdgePtr> edges;
        for (const auto& f : faces_) {
            for (const auto& e : *f) {
                edges.emplace_back(e);
            }
        }
        return edges;
    }

    /** @brief Get the list of faces in insertion order */
    auto faces() const -> std::vector<FacePtr> { return faces_; }

    /** @brief Get the list of interior vertices in insertion order */
    auto vertices_interior() const -> std::vector<VertPtr>
    {
        std::vector<VertPtr> ret;
        std::copy_if(
            verts_.begin(), verts_.end(), std::back_inserter(ret),
            [](auto x) { return not x->is_boundary(); });
        return ret;
    }

    /** @brief Get the list of boundary vertices in insertion order */
    auto vertices_boundary() const -> std::vector<VertPtr>
    {
        std::vector<VertPtr> ret;
        std::copy_if(
            verts_.begin(), verts_.end(), std::back_inserter(ret),
            [](auto x) { return x->is_boundary(); });
        return ret;
    }

    /** @brief Get the number of vertices */
    auto num_vertices() const -> std::size_t { return verts_.size(); }

    /** @brief Get the number of interior vertices */
    auto num_vertices_interior() const -> std::size_t
    {
        return std::accumulate(
            verts_.begin(), verts_.end(), std::size_t{0}, [](auto a, auto b) {
                return a + static_cast<std::size_t>(not b->is_boundary());
            });
    }

    /** @brief Get the number of edges */
    auto num_edges() const -> std::size_t { return edges_.size(); }

    /** @brief Get the number of faces */
    auto num_faces() const -> std::size_t { return faces_.size(); }

private:
    /** Find an existing edge with the provided end points */
    auto find_edge_(std::size_t start, std::size_t end) -> EdgePtr
    {
        // Get edges with this start index
        auto range = edges_.equal_range(start);

        // Loop over potential edges
        for (auto it = range.first; it != range.second; it++) {
            const auto& e = it->second;
            if (e->next and e->next->vertex->idx == end) {
                return e;
            }
        }
        return nullptr;
    }
};
}  // namespace OpenABF