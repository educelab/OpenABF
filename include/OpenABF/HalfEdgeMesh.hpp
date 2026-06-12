#pragma once

#include <algorithm>
#include <array>
#include <iterator>
#include <map>
#include <memory>
#include <optional>
#include <queue>
#include <sstream>
#include <unordered_set>
#include <vector>

#include "OpenABF/Exceptions.hpp"
#include "OpenABF/Vec.hpp"

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

namespace detail
{
/** Debug: Print a vector of elements to a string */
template <typename T>
auto vec_to_string(const T& v) -> std::string
{
    std::ostringstream ss;
    ss << '[';
    for (std::size_t i = 0; i < std::size(v); ++i) {
        if (i != 0)
            ss << ", ";
        ss << std::begin(v)[i];
    }
    ss << ']';
    return ss.str();
}

/**
 * Returns a copy of @p v with every element matching @p p removed. Naming
 * matches std::remove_if (predicate selects elements to drop), NOT the
 * conventional `filter` semantic where the predicate selects what to keep.
 */
template <class ForwardContainer, class UnaryPred>
auto remove_if(ForwardContainer v, UnaryPred p)
{
    auto end = std::remove_if(std::begin(v), std::end(v), p);
    v.erase(end, std::end(v));
    return v;
}

/**
 * @brief Lightweight non-owning range adapter holding begin/end iterators
 *
 * Supports range-based for loops, empty(), and front(). Used as the return
 * type for lazy mesh iteration methods.
 */
template <typename Iter>
struct Range {
    Iter first_;
    Iter last_;
    auto begin() const -> Iter { return first_; }
    auto end() const -> Iter { return last_; }
    auto empty() const -> bool { return first_ == last_; }
    auto front() const -> decltype(*first_) { return *first_; }
};

/**
 * @brief Input iterator that skips elements not matching a predicate
 *
 * Wraps a base iterator and advances automatically past elements for which
 * pred(*it) is false.
 */
template <typename Iter, typename Pred>
class FilteringIterator
{
public:
    using difference_type = std::ptrdiff_t;
    using value_type = typename std::iterator_traits<Iter>::value_type;
    using pointer = typename std::iterator_traits<Iter>::pointer;
    using reference = typename std::iterator_traits<Iter>::reference;
    using iterator_category = std::input_iterator_tag;

    FilteringIterator() = default;
    FilteringIterator(Iter current, Iter end, Pred pred) : current_{current}, end_{end}, pred_{pred}
    {
        advance_to_next();
    }

    auto operator*() const -> reference { return *current_; }
    auto operator->() const -> pointer { return &*current_; }

    auto operator++() -> FilteringIterator&
    {
        ++current_;
        advance_to_next();
        return *this;
    }

    auto operator==(const FilteringIterator& other) const -> bool
    {
        return current_ == other.current_;
    }
    auto operator!=(const FilteringIterator& other) const -> bool { return !(*this == other); }

private:
    void advance_to_next()
    {
        while (current_ != end_ && !pred_(*current_)) {
            ++current_;
        }
    }

    Iter current_{};
    Iter end_{};
    Pred pred_{};
};
}  // namespace detail

/**
 * @brief Compute the internal angles of a face
 *
 * Updates the current angle (DefaultEdgeTraits::alpha) with the internal angles
 * derived from the face's vertex positions. Useful if you want to reset a face
 * after being processed by ABF or ABFPlusPlus.
 *
 * @tparam FacePtr A Face-type pointer implementing DefaultEdgeTraits
 * @throws MeshException If interior angle is NaN or Inf
 */
template <class FacePtr>
void ComputeFaceAngles(FacePtr& face)
{
    for (auto& e : *face) {
        auto ab = e->next->vertex->pos - e->vertex->pos;
        auto ac = e->next->next->vertex->pos - e->vertex->pos;
        e->alpha = interior_angle(ab, ac);
        if (std::isnan(e->alpha) or std::isinf(e->alpha)) {
            auto msg = "Interior angle for edge " + std::to_string(e->idx) + " is nan/inf";
            throw MeshException(msg);
        }
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

/** @brief Determines if mesh is open or closed */
template <class MeshPtr>
auto HasBoundary(const MeshPtr& mesh) -> bool
{
    for (const auto& v : mesh->vertices()) {
        if (v->is_boundary()) {
            return true;
        }
    }
    return false;
}

/**
 * @brief Check if a mesh has unreferenced vertices
 *
 * @note This only checks if the vertex is associated with at least one edge.
 * A face is not currently guaranteed.
 */
template <class MeshPtr>
auto HasUnreferencedVertices(const MeshPtr& mesh) -> bool
{
    for (const auto& v : mesh->vertices()) {
        if (v->edges.empty()) {
            return true;
        }
    }
    return false;
}

/**
 * @brief Get a list of unreferenced vertices
 *
 * @note This only checks if the vertex is associated with at least one edge.
 * A face is not currently guaranteed.
 */
template <class MeshPtr>
auto UnreferencedVertices(const MeshPtr& mesh) -> std::vector<std::size_t>
{
    std::vector<std::size_t> indices;
    for (const auto& v : mesh->vertices()) {
        if (v->edges.empty()) {
            indices.emplace_back(v->idx);
        }
    }
    return indices;
}

/** @brief Check if mesh is manifold */
template <class MeshPtr>
auto IsManifold(const MeshPtr& mesh) -> bool
{
    // insert_faces won't allow non-manifold edges, but vertices may still be
    // non-manifold if update_boundary was never called, so check those here
    for (const auto& v : mesh->vertices()) {
        if (not v->is_manifold()) {
            return false;
        }
    }
    return true;
}

/**
 * @brief Find an edge path between two vertices
 *
 * Uses Dijkstra's algorithm to find the shortest path between two vertices.
 * Distance is measured using the edge lengths of the mesh. The returned mesh
 * is not guaranteed to be the _only_ shortest path (there may be many which
 * have the same length), but only the first discovered.
 *
 * If the returned list is empty, the endpoints are the same or a path between
 * the two endpoints does not exist (i.e. the mesh has multiple connected
 * components).
 *
 * @returns std::vector<EdgePtr>
 */
template <class MeshPtr>
auto FindEdgePath(const MeshPtr& mesh, std::size_t from, std::size_t to)
{
    using Mesh = std::remove_reference_t<decltype(*mesh)>;
    using Value = typename Mesh::type;
    using EdgePtr = typename Mesh::EdgePtr;

    // End points are the same
    if (from == to) {
        return std::vector<EdgePtr>{};
    }

    struct Node {
        using Ptr = std::shared_ptr<Node>;
        std::size_t idx{0};
        Value dist{INF<Value>};
        Ptr prev{nullptr};
        EdgePtr fromEdge{nullptr};
    };

    // Build a list of all nodes
    std::vector<typename Node::Ptr> nodes;
    nodes.reserve(mesh->num_vertices());
    for (std::size_t i = 0; i < mesh->num_vertices(); ++i) {
        auto n = std::make_shared<Node>();
        n->idx = i;
        if (i == from) {
            n->dist = 0;
        }
        nodes.push_back(n);
    }

    // Build a queue
    struct Compare {
        auto operator()(const typename Node::Ptr& p, const typename Node::Ptr& q) const -> bool
        {
            return p->dist > q->dist;
        }
    };
    using Queue = std::priority_queue<typename Node::Ptr, std::vector<typename Node::Ptr>, Compare>;
    Queue queue;
    queue.push(nodes[from]);

    typename Node::Ptr end{nullptr};
    while (not queue.empty()) {
        auto p = queue.top();
        queue.pop();
        if (p->idx == to) {
            end = p;
            break;
        }

        for (const auto& e : mesh->outgoing_edges(p->idx)) {
            const auto next = e->pair->vertex->idx;
            const auto d = p->dist + e->magnitude();

            if (d < nodes[next]->dist) {
                nodes[next]->dist = d;
                nodes[next]->prev = p;
                nodes[next]->fromEdge = e;
                queue.push(nodes[next]);
            }
        }
    }

    // Error: haven't found a path
    if (end == nullptr) {
        return std::vector<EdgePtr>{};
    }

    // Build edge path
    std::vector<EdgePtr> path;
    auto node = end;
    while (node->prev) {
        path.emplace_back(node->fromEdge);
        node = node->prev;
    }
    std::reverse(path.begin(), path.end());
    return path;
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
template <typename T, std::size_t Dim = 3, typename VertexTraits = traits::DefaultVertexTraits<T>,
          typename EdgeTraits = traits::DefaultEdgeTraits<T>,
          typename FaceTraits = traits::DefaultFaceTraits<T>>
class HalfEdgeMesh
{
public:
    /** Fundamental element type (e.g. float, double) */
    using type = T;

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
        using pointer = std::conditional_t<Const, value_type const*, value_type*>;
        /** Reference type */
        using reference = std::conditional_t<Const, value_type const&, value_type&>;
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
        auto operator!=(const FaceIterator& other) const -> bool { return !(*this == other); }
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

    /**
     * @brief Iterator for the edges of a vertex's wheel
     *
     * Walks the half-edge linked-list around a vertex (`edge->pair->next`),
     * yielding only non-boundary edges. Multi-pass safe: constructing a new
     * WheelIterator from the same vertex edge always starts at the beginning.
     *
     * @tparam Const If true, is a const iterator
     */
    template <bool Const = false>
    class WheelIterator
    {
    public:
        /** Difference type */
        using difference_type = std::ptrdiff_t;
        /** Value type */
        using value_type = EdgePtr;
        /** Pointer type */
        using pointer = std::conditional_t<Const, value_type const*, value_type*>;
        /** Reference type */
        using reference = std::conditional_t<Const, value_type const&, value_type&>;
        /** Iterator category */
        using iterator_category = std::input_iterator_tag;

        /** Default constructor == End iterator (current_ == nullptr) */
        WheelIterator() = default;
        /** Construct from the vertex's stored edge */
        explicit WheelIterator(const EdgePtr& head) : head_{head}, current_{head}
        {
            advance_to_non_boundary();
        }

        /** Dereference */
        template <bool C = Const>
        auto operator*() const -> std::enable_if_t<C, reference>
        {
            return current_;
        }
        template <bool C = Const>
        auto operator*() -> std::enable_if_t<!C, reference>
        {
            return current_;
        }

        /** Equality */
        auto operator==(const WheelIterator& other) const -> bool
        {
            return current_ == other.current_;
        }
        /** Inequality */
        auto operator!=(const WheelIterator& other) const -> bool { return !(*this == other); }

        /** Increment */
        auto operator++() -> WheelIterator&
        {
            if (!current_) {
                return *this;
            }
            current_ = current_->pair->next;
            if (current_ == head_) {
                current_ = nullptr;
                return *this;
            }
            advance_to_non_boundary();
            return *this;
        }

    private:
        void advance_to_non_boundary()
        {
            while (current_ && current_->is_boundary()) {
                current_ = current_->pair->next;
                if (current_ == head_) {
                    current_ = nullptr;
                    break;
                }
            }
        }

        EdgePtr head_{};
        EdgePtr current_{};
    };

    /**
     * @brief Iterator that lazily flattens all face edges across all faces
     *
     * Provides a single flat sequence over all face (non-boundary) edges in
     * the mesh without allocating a vector. Internally advances through
     * `faces_` and uses `FaceIterator` within each face.
     *
     * @tparam Const If true, is a const iterator
     */
    template <bool Const = false>
    class EdgesIterator
    {
    public:
        /** Difference type */
        using difference_type = std::ptrdiff_t;
        /** Value type */
        using value_type = EdgePtr;
        /** Pointer type */
        using pointer = std::conditional_t<Const, value_type const*, value_type*>;
        /** Reference type */
        using reference = std::conditional_t<Const, value_type const&, value_type&>;
        /** Iterator category */
        using iterator_category = std::input_iterator_tag;

        using FaceVecIter = typename std::vector<FacePtr>::const_iterator;

        /** Construct at position (begin or end depending on faceIt == faceEnd) */
        EdgesIterator(FaceVecIter faceIt, FaceVecIter faceEnd) : faceIt_{faceIt}, faceEnd_{faceEnd}
        {
            if (faceIt_ != faceEnd_) {
                edgeIt_ = FaceIterator<true>{(*faceIt_)->head, (*faceIt_)->head};
                advance_if_face_exhausted();
            }
        }

        /** Dereference */
        template <bool C = Const>
        auto operator*() const -> std::enable_if_t<C, reference>
        {
            return *edgeIt_;
        }
        template <bool C = Const>
        auto operator*() -> std::enable_if_t<!C, reference>
        {
            return *edgeIt_;
        }

        /** Equality */
        auto operator==(const EdgesIterator& other) const -> bool
        {
            // Both exhausted (at end)
            if (faceIt_ == faceEnd_ && other.faceIt_ == other.faceEnd_) {
                return true;
            }
            if (faceIt_ != other.faceIt_) {
                return false;
            }
            return edgeIt_ == other.edgeIt_;
        }
        /** Inequality */
        auto operator!=(const EdgesIterator& other) const -> bool { return !(*this == other); }

        /** Increment */
        auto operator++() -> EdgesIterator&
        {
            ++edgeIt_;
            advance_if_face_exhausted();
            return *this;
        }

    private:
        void advance_if_face_exhausted()
        {
            while (edgeIt_ == FaceIterator<true>() && faceIt_ != faceEnd_) {
                ++faceIt_;
                if (faceIt_ != faceEnd_) {
                    edgeIt_ = FaceIterator<true>{(*faceIt_)->head, (*faceIt_)->head};
                }
            }
        }

        FaceVecIter faceIt_{};
        FaceVecIter faceEnd_{};
        FaceIterator<true> edgeIt_{};
    };

public:
    /** @brief %Vertex type */
    struct Vertex : VertexTraits {
        /** @brief Default constructor */
        Vertex() = default;

        /** @brief Construct from position values */
        template <typename... Args>
        explicit Vertex(Args... args) : pos{args...}
        {
        }

        /** @brief Copy-construct inherited traits */
        Vertex(const Vertex& rhs) : VertexTraits(rhs), pos{rhs.pos} {}

        /** @brief Construct a new Vertex pointer */
        template <typename... Args>
        static auto New(Args&&... args) -> VertPtr
        {
            return std::make_shared<Vertex>(std::forward<Args>(args)...);
        }

        /**
         * @brief Get the edges of a vertex's wheel
         *
         * @throws MeshException If vertex is a boundary vertex.
         */
        auto wheel() const
        {
            using Iter = WheelIterator<true>;
            return detail::Range<Iter>{Iter{edge}, Iter{}};
        }

        /** @brief Unit vertex normal */
        auto normal() const -> Vec<T, Dim>
        {
            Vec<T, Dim> n{0, 0, 0};
            for (const auto& e : wheel()) {
                n += e->face->normal();
            }
            return normalize(n);
        }

        /** @brief Returns if vertex is on mesh boundary */
        [[nodiscard]] auto is_boundary() const -> bool
        {
            auto e = edge;
            do {
                if (e->is_boundary() or e->pair->is_boundary()) {
                    return true;
                }
                e = e->pair->next;
            } while (e != edge);
            return false;
        }

        /** @brief Returns if vertex is interior to mesh */
        [[nodiscard]] auto is_interior() const -> bool { return not is_boundary(); }

        /** @brief Returns if vertex is unreferenced */
        [[nodiscard]] auto is_unreferenced() const -> bool { return edge == nullptr; }

        /** @brief Returns if vertex is manifold */
        [[nodiscard]] auto is_manifold() const -> bool
        {
            std::size_t boundaryCnt{0};
            auto out = mesh->outgoing_edges(idx);
            for (const auto& e : out) {
                if (e->is_boundary() or e->pair->is_boundary()) {
                    boundaryCnt++;
                }
            }
            return boundaryCnt == 0 or boundaryCnt == 2;
        }

        /** @brief Insertion index */
        std::size_t idx{0};
        /** @brief Vertex position */
        Vec<T, Dim> pos;
        /**
         * @brief Pointer to _an_ Edge with this vertex as its head
         *
         * @note There may be many such vertices.
         */
        EdgePtr edge;
        /** @brief Mesh to which this vertex belongs */
        HalfEdgeMesh* mesh{nullptr};
    };

    /** @brief %Edge type */
    struct Edge : EdgeTraits {
        Edge() = default;

        /** @brief Copy-construct inherited traits */
        Edge(const Edge& rhs) : EdgeTraits(rhs) {}

        /** @brief Construct a new Edge pointer */
        template <typename... Args>
        static auto New(Args&&... args) -> EdgePtr
        {
            return std::make_shared<Edge>(std::forward<Args>(args)...);
        }

        /** @brief Returns if edge is on mesh boundary */
        [[nodiscard]] auto is_boundary() const -> bool { return face == nullptr; }

        /** @brief Edge length */
        auto magnitude() const -> T { return (pair->vertex->pos - vertex->pos).magnitude(); }

        /** @brief This edge's adjacent half-edge */
        EdgePtr pair;
        /**
         * @brief The next edge in this edge's face
         *
         * If the edge is not assigned to a face, the next edge along the
         * boundary.
         */
        EdgePtr next;
        /**
         * @brief The previous edge in this edge's face
         *
         * If the edge is not assigned to a face, the previous edge along the
         * boundary.
         */
        EdgePtr prev;
        /** @brief The edge's vertex */
        VertPtr vertex;
        /** @brief The face containing this edge */
        FacePtr face;
        /** @brief Insertion index among all edges (including boundary) */
        std::size_t idx{0};
        /** @brief Insertion index among all face edges (excluding boundary) */
        std::optional<std::size_t> idxI;
        /** @brief Mesh to which this edge belongs */
        HalfEdgeMesh* mesh{nullptr};
    };

    /** @brief %Face type */
    struct Face : FaceTraits {
        /** Default constructor */
        Face() = default;

        /** @brief Copy-construct inherited traits */
        Face(const Face& rhs) : FaceTraits(rhs) {}

        /** @brief Construct a new Face pointer */
        template <typename... Args>
        static auto New(Args&&... args) -> FacePtr
        {
            return std::make_shared<Face>(std::forward<Args>(args)...);
        }

        /** Face edge iterator type */
        using iterator = FaceIterator<>;
        /** Face edge const iterator type */
        using const_iterator = FaceIterator<true>;
        /** @brief Returns an iterator over the edges of the face */
        iterator begin() { return iterator{head, head}; }
        /** @brief Returns the end iterator */
        iterator end() { return iterator(); }
        /** @brief Returns a const iterator over the edges of the face */
        const_iterator cbegin() const { return const_iterator{head, head}; }
        /** @brief Returns the const end iterator */
        const_iterator cend() const { return const_iterator(); }
        /** @brief Returns an iterable range over the edges of the face */
        auto edges() const
        {
            using Iter = const_iterator;
            return detail::Range<Iter>{Iter{head, head}, Iter{}};
        }

        /** @brief Area of the face */
        auto area() const -> T
        {
            // Get the edge lengths
            std::array<T, 3> l{head->magnitude(), head->next->magnitude(), head->prev->magnitude()};

            // Sort the side lengths so that a >= b >= c
            std::sort(l.begin(), l.end(), std::greater<T>());

            // Calculate the area
            const auto& a = l[0];
            const auto& b = l[1];
            const auto& c = l[2];
            auto p = (a + (b + c)) * (c - (a - b)) * (c + (a - b)) * (a + (b - c));
            return 0.25 * std::sqrt(p);
        }

        /** @brief Face barycenter (center-of-mass) */
        auto barycenter() const -> Vec<T, Dim>
        {
            return (head->vertex->pos + head->next->vertex->pos + head->prev->vertex->pos) / T(3);
        }

        /** @brief Unit face normal */
        auto normal() const -> Vec<T, Dim>
        {
            // Get the edge vectors
            auto e0 = head->prev->vertex->pos - head->vertex->pos;
            auto e1 = head->next->vertex->pos - head->vertex->pos;

            // Take the cross-product
            return normalize(e1.cross(e0));
        }

        /** @brief First edge in the face */
        EdgePtr head;
        /** @brief The next face in the mesh */
        FacePtr next;
        /** @brief Insertion index */
        std::size_t idx{0};
        /** @brief Mesh to which this vertex belongs */
        HalfEdgeMesh* mesh{nullptr};
    };

private:
    /** List of vertices */
    std::vector<VertPtr> verts_;
    /** List of faces */
    std::vector<FacePtr> faces_;
    /** List of all edges, indexed by the vertex's insertion index */
    std::multimap<std::size_t, EdgePtr> edges_;
    /** Number of all edges which border a face */
    std::size_t numFaceEdges_{0};

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
            e.second->prev = nullptr;
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
     * @brief Clone this mesh
     *
     * Returns a new mesh with the same structure as this mesh but not sharing
     * vertex, face, or edge elements.
     */
    auto clone() const -> Pointer
    {
        auto ret = HalfEdgeMesh::New();
        ret->verts_.reserve(verts_.size());
        ret->faces_.reserve(faces_.size());

        // Insert vertices
        for (const auto& v : verts_) {
            auto i = ret->insert_vertex(*v);
            ret->verts_[i]->edge = nullptr;
        }

        // Insert faces and edges
        for (const auto& f : faces_) {
            ret->clone_face_(f);
        }
        ret->update_boundary();
        return ret;
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
        vert->mesh = this;
        vert->idx = verts_.size();
        verts_.push_back(vert);
        return vert->idx;
    }

    /**
     * @brief Insert new vertices from a list of Vertex-like objects
     *
     * A convenience function which adds multiple vertices to the mesh.
     */
    template <class VectorOfVectors>
    auto insert_vertices(const VectorOfVectors& v) -> std::vector<std::size_t>
    {
        std::vector<std::size_t> idxs;
        for (const auto& f : v) {
            idxs.emplace_back(insert_vertex(f));
        }
        return idxs;
    }

    /**
     * @copydoc insert_vertices(const VectorOfVectors&)
     */
    template <typename ValType>
    auto insert_vertices(std::initializer_list<std::initializer_list<ValType>> v)
        -> std::vector<std::size_t>
    {
        auto it = std::begin(v);
        std::vector<std::size_t> idxs;
        for (std::size_t i = 0; i < std::size(v); ++i) {
            idxs.emplace_back(insert_vertex(it[i]));
        }
        return idxs;
    }

    /**
     * @brief Insert a face from an ordered list of Vertex indices
     *
     * Accepts an iterable supporting range-based for loops.
     *
     * @note This function does **not** update the mesh boundary connections.
     * Call update_boundary() after all faces have been inserted, or use
     * insert_faces() to insert faces and update the boundary in one step.
     *
     * @param vector List of vertex indices
     * @throws std::out_of_range If one of the vertex indices is out of bounds.
     * @throws MeshException (1) If one of provided edges is already paired.
     * This indicates that the mesh is not 2-manifold. (2) If an edge has
     * zero length. This means the face has zero area. (3) If an edge's interior
     * angle is NaN or Inf.
     */
    template <class Vector>
    auto insert_face(Vector&& vector) -> std::size_t
    {
        return insert_face_(std::forward<Vector>(vector));
    }

    /**
     * @brief Insert a new face from an ordered list of Vertex indices
     *
     * Accepts vertex indices as individual variadic arguments.
     *
     * @note This function does **not** update the mesh boundary connections.
     * Call update_boundary() after all faces have been inserted, or use
     * insert_faces() to insert faces and update the boundary in one step.
     *
     * @throws std::out_of_range If one of the vertex indices is out of bounds.
     * @throws MeshException If one of provided edges is already paired. This
     * indicates that the mesh is not 2-manifold.
     */
    template <typename... Args>
    auto insert_face(Args... args) -> std::size_t
    {
        static_assert(sizeof...(args) >= 3, "Faces require >= 3 indices");
        using Tuple = std::tuple<Args...>;
        using ElemT = std::tuple_element_t<0, Tuple>;
        return insert_face_(std::initializer_list<ElemT>{args...});
    }

    /**
     * @brief Insert new faces from a list of lists of Vertex indices
     *
     * A convenience function which adds multiple faces to the mesh (using
     * insert_face()) and updates the mesh boundary (using update_boundary())
     * when complete.
     */
    template <class VectorOfVectors>
    auto insert_faces(const VectorOfVectors& v) -> std::vector<std::size_t>
    {
        std::vector<std::size_t> idxs;
        for (const auto& f : v) {
            idxs.emplace_back(insert_face_(f));
        }

        update_boundary();
        return idxs;
    }

    /**
     * @copydoc insert_faces()
     */
    template <typename IdxType>
    auto insert_faces(std::initializer_list<std::initializer_list<IdxType>> v)
    {
        const auto it = std::begin(v);
        std::vector<std::size_t> idxs;
        for (std::size_t i = 0; i < std::size(v); ++i) {
            idxs.emplace_back(insert_face(it[i]));
        }
        update_boundary();
        return idxs;
    }

    /**
     * @brief Update the mesh boundary connections
     *
     * Because the mesh boundary may become temporarily non-traversable while
     * the mesh is being constructed, the mesh boundary connections should only
     * be updated after all faces have been added to the mesh. Call this
     * function after inserting faces with insert_face() or use insert_faces()
     * to construct the mesh and update the boundary in one step.
     */
    void update_boundary()
    {
        for (const auto& [_, edge] : edges_) {
            if (edge->is_boundary()) {
                // Get incoming boundary edges to the start point
                auto inBoundary =
                    detail::remove_if(incoming_edges(edge->vertex->idx),
                                      [](const auto& e) { return not e->is_boundary(); });
                if (inBoundary.size() == 0 or inBoundary.size() > 1) {
                    const std::array<std::size_t, 2> idx{edge->vertex->idx,
                                                         edge->pair->vertex->idx};
                    throw MeshException("Cannot update mesh boundary along edge " +
                                        detail::vec_to_string(idx) +
                                        " due to non-manifold surface and/or inconsistent "
                                        "winding order");
                }

                // Get outgoing boundary edges to the end point
                auto outBoundary =
                    detail::remove_if(outgoing_edges(edge->pair->vertex->idx),
                                      [](const auto& e) { return not e->is_boundary(); });
                if (outBoundary.size() == 0 or outBoundary.size() > 1) {
                    const std::array<std::size_t, 2> idx{edge->vertex->idx,
                                                         edge->pair->vertex->idx};
                    throw MeshException("Cannot update mesh boundary along edge " +
                                        detail::vec_to_string(idx) +
                                        " due to non-manifold surface and/or inconsistent "
                                        "winding order");
                }

                edge->prev = inBoundary[0];
                inBoundary[0]->next = edge;
                edge->next = outBoundary[0];
                outBoundary[0]->prev = edge;
            }
        }
    }

    /** @brief Get the list of vertices in insertion order */
    auto vertices() const -> const std::vector<VertPtr>& { return verts_; }

    /** @brief Get a vertex by index */
    auto vertex(std::size_t idx) const -> VertPtr { return verts_.at(idx); }

    /** @brief Get a lazy range over all face edges in insertion order */
    auto edges() const
    {
        using Iter = EdgesIterator<true>;
        return detail::Range<Iter>{Iter{faces_.cbegin(), faces_.cend()},
                                   Iter{faces_.cend(), faces_.cend()}};
    }

    /** @brief Find an existing edge with the provided end points */
    auto edge(std::size_t start, std::size_t end) -> EdgePtr
    {
        // Get edges with this start index
        const auto range = edges_.equal_range(start);

        // Loop over potential edges
        for (auto it = range.first; it != range.second; ++it) {
            const auto& e = it->second;
            if (e->pair->vertex->idx == end) {
                return e;
            }
        }
        return nullptr;
    }

    /**
     * @brief Get a boundary edge
     *
     * Returns the first boundary edge in the list of edges
     */
    auto boundary_edge() const -> EdgePtr
    {
        // Find a boundary edge
        for (const auto& e : edges_) {
            if (e.second->is_boundary()) {
                return e.second;
            }
        }
        return nullptr;
    }

    /**
     * @brief Build a list of all boundary edges
     *
     * Returns a list of lists of edges which lie on one of this mesh's
     * boundaries. One edge list is returned for each unique boundary. Meshes
     * with a single connected component may still have multiple boundaries
     * (i.e. if the mesh has holes).
     */
    auto boundaries() const -> std::vector<std::vector<EdgePtr>>
    {
        using Boundary = std::vector<EdgePtr>;
        Boundary boundary;
        std::vector<Boundary> boundaries;
        std::unordered_set<std::size_t> visited;
        for (const auto& e : edges_) {
            const auto edge = e.second;
            // Skip is not a boundary edge or already visited
            if (not edge->is_boundary() or visited.count(edge->idx) > 0) {
                continue;
            }

            // Create a new boundary
            boundary.clear();
            visited.insert(edge->idx);
            boundary.emplace_back(edge);
            auto test = edge->next;
            do {
                visited.insert(test->idx);
                boundary.emplace_back(test);
                test = test->next;
            } while (test != edge);
            boundaries.emplace_back(boundary);
        }
        return boundaries;
    }

    /** @brief Get the list of faces in insertion order */
    auto faces() const -> const std::vector<FacePtr>& { return faces_; }

    /** @brief Get a face by index */
    auto face(std::size_t idx) const -> FacePtr { return faces_.at(idx); }

    /**
     * @brief Get the number of connected components
     *
     * A connected component is a set of continuous, adjacent faces. If you
     * want to get the list of connected components, use connected_components().
     *
     * @see connected_components()
     */
    [[nodiscard]] auto num_connected_components() const -> std::size_t
    {
        std::size_t cnt{0};
        std::vector visited(num_faces(), false);
        std::queue<FacePtr> queue;
        // Iterate over the faces
        for (const auto& f : faces_) {
            // Skip faces we've visited
            if (visited[f->idx]) {
                continue;
            }

            // Start a new connected component
            queue.push(f);
            while (not queue.empty()) {
                // Get the top of the queue
                auto p = queue.front();
                queue.pop();
                // Mark as visited
                visited[p->idx] = true;
                // Add the neighbor faces to the queue
                for (const auto& e : *p) {
                    if (not e->pair->is_boundary()) {
                        auto n = e->pair->face;
                        if (not visited[n->idx]) {
                            queue.push(n);
                        }
                    }
                }
            }
            // Finished this component
            ++cnt;
        }
        return cnt;
    }

    /** @brief Get a list of connected components */
    auto connected_components() const -> std::vector<std::vector<FacePtr>>
    {
        // Tracking structures
        std::vector<std::vector<FacePtr>> components;
        std::vector visited(num_faces(), false);
        std::vector<FacePtr> current;
        std::queue<FacePtr> queue;

        // Iterate over the faces
        for (const auto& f : faces_) {
            // Skip faces we've visited
            if (visited[f->idx]) {
                continue;
            }

            // Start a new connected component
            current.clear();
            queue.push(f);
            while (not queue.empty()) {
                // Get the top of the queue
                auto p = queue.front();
                queue.pop();
                // Mark as visited
                visited[p->idx] = true;
                // Add to this connected component
                current.emplace_back(p);
                // Add the neighbor faces to the queue
                for (const auto& e : *p) {
                    if (not e->pair->is_boundary()) {
                        auto n = e->pair->face;
                        if (not visited[n->idx]) {
                            queue.push(n);
                        }
                    }
                }
            }
            // Add this component to the list
            components.emplace_back(current);
        }
        return components;
    }

    /** @brief Get the list of interior vertices in insertion order */
    auto vertices_interior() const
    {
        auto pred = [](const VertPtr& v) { return not v->is_boundary(); };
        using BaseIter = typename std::vector<VertPtr>::const_iterator;
        using Iter = detail::FilteringIterator<BaseIter, decltype(pred)>;
        return detail::Range<Iter>{Iter{verts_.cbegin(), verts_.cend(), pred},
                                   Iter{verts_.cend(), verts_.cend(), pred}};
    }

    /** @brief Get a lazy range over boundary vertices in insertion order */
    auto vertices_boundary() const
    {
        auto pred = [](const VertPtr& v) { return v->is_boundary(); };
        using BaseIter = typename std::vector<VertPtr>::const_iterator;
        using Iter = detail::FilteringIterator<BaseIter, decltype(pred)>;
        return detail::Range<Iter>{Iter{verts_.cbegin(), verts_.cend(), pred},
                                   Iter{verts_.cend(), verts_.cend(), pred}};
    }

    /** @brief Get the number of vertices */
    [[nodiscard]] auto num_vertices() const -> std::size_t { return verts_.size(); }

    /** @brief Get the number of interior vertices */
    [[nodiscard]] auto num_vertices_interior() const -> std::size_t
    {
        return std::accumulate(verts_.begin(), verts_.end(), std::size_t{0}, [](auto a, auto b) {
            return a + static_cast<std::size_t>(not b->is_boundary());
        });
    }

    /** @brief Get the number of edges */
    [[nodiscard]] auto num_edges() const -> std::size_t
    {
        std::size_t ret = 0;
        for (const auto& [_, e] : edges_) {
            if (not e->is_boundary()) {
                ++ret;
            }
        }
        return ret;
    }

    /** @brief Get the number of faces */
    [[nodiscard]] auto num_faces() const -> std::size_t { return faces_.size(); }

    /**
     * @brief Split an edge in order to introduce a new boundary
     *
     * Disconnects a single, paired half-edge (i.e. a single edge between two
     * connected triangles) into two, unpaired half-edges, creating a
     * "hole" in the mesh's connectivity graph. Useful when you want to
     * introduce a boundary, or tear, into a mesh to improve parameterization.
     *
     * @note If the given edge intersects with an existing boundary, one
     * or both of your endpoint vertices will be duplicated. Be sure to take
     * this into account when converting the parameterized mesh to your
     * (per-wedge) UV map.
     *
     * @see split_path()
     */
    void split_edge(const EdgePtr& edge)
    {
        // Get forward and backward edge
        auto oldFwd = edge;
        auto oldBwd = oldFwd->pair;

        // Don't split boundary edge pairs
        if (oldFwd->is_boundary() and oldBwd->is_boundary()) {
            return;
        }

        // Get initial vertices
        auto oldStart = oldFwd->vertex;
        auto oldEnd = oldBwd->vertex;
        auto startOnBoundary = oldStart->is_boundary();
        auto endOnBoundary = oldEnd->is_boundary();

        // Get the new starting vertex for this edge. detail::remove_if keeps
        // only boundary half-edges (predicate selects what to drop). A
        // manifold boundary vertex has exactly one boundary half-edge in each
        // direction; more than one means a non-manifold split endpoint.
        VertPtr newStart;
        EdgePtr startIn, startOut;
        if (startOnBoundary) {
            auto newIdx = insert_vertex(oldStart->pos);
            newStart = verts_.at(newIdx);

            auto in = detail::remove_if(incoming_edges(oldStart->idx),
                                        [](auto e) { return not e->is_boundary(); });
            auto out = detail::remove_if(outgoing_edges(oldStart->idx),
                                         [](auto e) { return not e->is_boundary(); });
            if (in.size() != 1 or out.size() != 1) {
                throw MeshException("Non-manifold boundary vertex at split endpoint");
            }
            startIn = in[0];
            startOut = out[0];
        } else {
            newStart = oldStart;
        }

        // Get the new ending vertex for this edge
        VertPtr newEnd;
        EdgePtr endIn, endOut;
        if (endOnBoundary) {
            auto newIdx = insert_vertex(oldEnd->pos);
            newEnd = verts_.at(newIdx);

            auto in = detail::remove_if(incoming_edges(oldEnd->idx),
                                        [](auto e) { return not e->is_boundary(); });
            auto out = detail::remove_if(outgoing_edges(oldEnd->idx),
                                         [](auto e) { return not e->is_boundary(); });
            if (in.size() != 1 or out.size() != 1) {
                throw MeshException("Non-manifold boundary vertex at split endpoint");
            }
            endIn = in[0];
            endOut = out[0];
        } else {
            newEnd = oldEnd;
        }

        // Create new edge pair
        auto newFwd = Edge::New();
        auto newBwd = Edge::New();
        newFwd->pair = newBwd;
        newBwd->pair = newFwd;
        newFwd->mesh = newBwd->mesh = this;

        // Assign vertices and add to mesh
        newFwd->vertex = newStart;
        newFwd->idx = edges_.size();
        edges_.emplace(newStart->idx, newFwd);
        newBwd->vertex = newEnd;
        newBwd->idx = edges_.size();
        edges_.emplace(newEnd->idx, newBwd);

        // Update vertices with edge if required
        if (not newStart->edge) {
            newStart->edge = newFwd;
        }
        if (not newEnd->edge) {
            newEnd->edge = newBwd;
        }

        // New forward takes old forward's face edge idx
        std::swap(newFwd->idxI, oldFwd->idxI);

        // Move old forward's face to new forward
        std::swap(newFwd->face, oldFwd->face);
        std::swap(newFwd->next, oldFwd->next);
        std::swap(newFwd->prev, oldFwd->prev);
        std::swap(newFwd->alpha, oldFwd->alpha);
        if (newFwd->face->head == oldFwd) {
            newFwd->face->head = newFwd;
        }

        // Re-link neighbors of newFwd in the inherited face.
        //
        // newFwd->next (second half-edge of the face) needs to originate at
        // newEnd, which it does unconditionally below — a no-op when
        // endOnBoundary is false (newEnd == oldEnd).
        //
        // The neighboring half-edge across newFwd->prev (newFwd->prev->pair)
        // also originates at newStart in the boundary branch, but it is
        // ALWAYS covered downstream: it is either startOut itself (when the
        // existing boundary at oldStart sits across newFwd->prev — the
        // SimpleSplit case) and updated by the `startOut->vertex = newStart`
        // assignment, or it is a non-boundary edge in newStart's fan and
        // updated by the wheel loop below.
        newFwd->next->prev = newFwd;
        newFwd->prev->next = newFwd;
        rekey_edge_to_vertex(newFwd->next, newEnd);

        // Update new boundary edges' next/prev. The rekey_* calls move
        // half-edges between edges_ multimap buckets so that
        // outgoing_edges(idx) / Vertex::is_manifold() remain consistent.
        if (startOnBoundary) {
            rekey_edge_to_vertex(startOut, newStart);
            newBwd->next = startOut;
            startOut->prev = newBwd;
            oldFwd->prev = startIn;
            startIn->next = oldFwd;
            for (auto e : newStart->wheel()) {
                rekey_edge_to_vertex(e, newStart);
            }
        } else {
            newBwd->next = oldFwd;
            oldFwd->prev = newBwd;
        }
        if (endOnBoundary) {
            newBwd->prev = endIn;
            endIn->next = newBwd;
            oldFwd->next = endOut;
            endOut->prev = oldFwd;
            for (auto e : newEnd->wheel()) {
                rekey_edge_to_vertex(e, newEnd);
            }
        } else {
            newBwd->prev = oldFwd;
            oldFwd->next = newBwd;
        }
    }

    /**
     * @brief Split a list of edges (path) to form a new boundary
     *
     * @note This function was designed to split a list of continuous edges
     * forming a path on the surface of the mesh. Providing otherwise can lead
     * to undefined behavior.
     *
     * @see split_edge()
     */
    void split_path(const std::vector<EdgePtr>& path)
    {
        // Split edges
        for (const auto& e : path) {
            split_edge(e);
        }
    }

    /**
     * @copybrief split_path
     *
     * This overload accepts a list of vertex indices forming a continuous path.
     *
     * @copydetails split_path
     */
    void split_path(const std::vector<std::size_t>& path)
    {
        // Convert index path to edge path
        std::vector<EdgePtr> edgePath;
        for (std::size_t i = 0; i < path.size() - 1; ++i) {
            auto e = this->edge(path[i], path[i + 1]);
            if (not e) {
                throw MeshException("Could not find edge");
            }
            edgePath.emplace_back(e);
        }

        split_path(edgePath);
    }

    /**
     * @brief Move @p e from its current edges_ bucket to the bucket keyed by
     * @p newVert and update e->vertex.
     *
     * The edges_ multimap is keyed by half-edge origin vertex index. Any code
     * that reassigns a half-edge's origin (e.g. split_edge duplicating a
     * vertex and moving a fan of half-edges to the new vertex) must call this
     * to keep the multimap consistent, otherwise outgoing_edges() /
     * incoming_edges() / Vertex::is_manifold() return stale results.
     */
    void rekey_edge_to_vertex(const EdgePtr& e, const VertPtr& newVert)
    {
        if (e->vertex == newVert) {
            return;
        }
        const auto range = edges_.equal_range(e->vertex->idx);
        for (auto it = range.first; it != range.second; ++it) {
            if (it->second == e) {
                edges_.erase(it);
                break;
            }
        }
        edges_.emplace(newVert->idx, e);
        e->vertex = newVert;
    }

    /** @brief Get a list of outgoing edges from a specific vertex (by index) */
    auto outgoing_edges(const std::size_t idx) -> std::vector<EdgePtr>
    {
        const auto range = edges_.equal_range(idx);
        std::vector<EdgePtr> ret;
        ret.reserve(std::distance(range.first, range.second));
        std::transform(range.first, range.second, std::back_inserter(ret),
                       [](auto it) { return it.second; });
        return ret;
    }

    /** @brief Get a list of incoming edges to a specific vertex (by index) */
    auto incoming_edges(const std::size_t idx) -> std::vector<EdgePtr>
    {
        auto outEdges = outgoing_edges(idx);
        std::vector<EdgePtr> ret;
        ret.reserve(outEdges.size());
        std::transform(outEdges.begin(), outEdges.end(), std::back_inserter(ret),
                       [](auto e) { return e->pair; });
        return ret;
    }

private:
    /**
     * Face insertion implementation
     *
     * @param vector Iterable of vertex indices
     * @param face Pre-existing Face (only used when cloning)
     */
    template <class Vector>
    auto insert_face_(const Vector& vector, FacePtr face = nullptr) -> std::size_t
    {
        // Make a new face structure
        if (not face) {
            face = Face::New();
        }
        face->mesh = this;

        // Create a list of vertex pairs
        using IDType = std::size_t;
        using IDPair = std::pair<std::size_t, IDType>;
        std::vector<IDPair> endPts;
        for (std::size_t i = 0; i < std::size(vector); ++i) {
            auto nextIdx = i == std::size(vector) - 1 ? 0 : i + 1;
            endPts.emplace_back(std::begin(vector)[i], std::begin(vector)[nextIdx]);
        }

        // Create a new edge for every edge pair
        bool reverse{false};
        std::vector<EdgePtr> edges;
        for (const auto& [startIdx, endIdx] : endPts) {
            // See if this edge already exists
            auto thisEdge = this->edge(startIdx, endIdx);

            // If this edge doesn't exist, make it and its pair
            if (not thisEdge) {
                thisEdge = Edge::New();
                auto pair = Edge::New();
                thisEdge->pair = pair;
                pair->pair = thisEdge;
                thisEdge->mesh = pair->mesh = this;

                thisEdge->idx = edges_.size();
                edges_.emplace(startIdx, thisEdge);
                pair->idx = edges_.size();
                edges_.emplace(endIdx, pair);

                thisEdge->vertex = verts_.at(startIdx);
                if (not thisEdge->vertex->edge) {
                    thisEdge->vertex->edge = thisEdge;
                }

                pair->vertex = verts_.at(endIdx);
                if (not pair->vertex->edge) {
                    pair->vertex->edge = pair;
                }
            }

            // Reverse winding order
            if (reverse) {
                thisEdge = thisEdge->pair;
            }

            // If this edge has a face, try reversing the winding order
            if (thisEdge->face) {
                auto pair = thisEdge->pair;
                // Winding order error if already reversed
                // TODO: Theoretically could recursively flip winding order for
                //       adjacent faces which violate the order
                if (reverse) {
                    const auto msg =
                        "Winding order cannot be fixed for face"
                        "with vids=" +
                        detail::vec_to_string(vector);
                    throw MeshException(msg);
                }
                // Non-manifold error if the pair is already assigned
                if (pair->face) {
                    const auto msg =
                        "Attempted to add non-manifold face along "
                        "edge with vids=" +
                        detail::vec_to_string(vector);
                    throw MeshException(msg);
                }
                reverse = true;
                thisEdge = pair;

                // If reversing, update existing visited edges with the pair
                for (auto& e : edges) {
                    pair = e->pair;
                    // If pair already has a face, then manifold error
                    if (pair->face) {
                        const auto msg =
                            "Attempted to add non-manifold face "
                            "along edge with vids=[" +
                            std::to_string(startIdx) + ", " + std::to_string(endIdx) + "]";
                        throw MeshException(msg);
                    }
                    e = pair;
                }
            }
            thisEdge->face = face;

            // Set the head edge for this face
            if (not face->head) {
                face->head = thisEdge;
            }

            // Store the edges for next/prev later
            edges.push_back(thisEdge);
        }

        // Update next/previous
        for (std::size_t i = 0; i < edges.size(); ++i) {
            auto edge = edges[i];
            const auto prevIdx = i == 0 ? edges.size() - 1 : i - 1;
            const auto nextIdx = i == edges.size() - 1 ? 0 : i + 1;
            edge->next = reverse ? edges[prevIdx] : edges[nextIdx];
            edge->prev = reverse ? edges[nextIdx] : edges[prevIdx];
            if (not edge->idxI) {
                edge->idxI = numFaceEdges_;
                ++numFaceEdges_;
            }
        }

        // Sanity check: edge lengths
        for (const auto& e : *face) {
            if (norm(e->next->vertex->pos - e->vertex->pos) == 0.0) {
                auto msg = "Zero-length edge (" + std::to_string(e->vertex->idx) + ", " +
                           std::to_string(e->next->vertex->idx) + ")";
                throw MeshException(msg);
            }
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

    /**
     * Extra steps which need to be run before insert_face_ when cloning a face
     * from an existing mesh
     *
     * @param face Existing face from the mesh being cloned
     */
    auto clone_face_(const FacePtr& face)
    {
        // Copy the existing face
        auto f = Face::New(*face);

        // Pre-make all edges
        std::vector<std::size_t> idxs;
        for (const auto& e : *face) {
            auto startIdx = e->vertex->idx;
            auto endIdx = e->pair->vertex->idx;
            idxs.emplace_back(startIdx);

            // Make sure we haven't created this edge and pair already
            auto outEdge = this->edge(startIdx, endIdx);

            // Copy all inherited properties
            if (not outEdge) {
                outEdge = Edge::New(*e);
                auto inEdge = Edge::New(*e->pair);

                outEdge->pair = inEdge;
                inEdge->pair = outEdge;
                outEdge->mesh = inEdge->mesh = this;

                outEdge->idx = edges_.size();
                edges_.emplace(startIdx, outEdge);
                inEdge->idx = edges_.size();
                edges_.emplace(endIdx, inEdge);

                outEdge->vertex = verts_.at(startIdx);
                if (not outEdge->vertex->edge) {
                    outEdge->vertex->edge = outEdge;
                }

                inEdge->vertex = verts_.at(endIdx);
                if (not inEdge->vertex->edge) {
                    inEdge->vertex->edge = inEdge;
                }
            }
        }

        // Create the face
        return insert_face_(idxs, f);
    }
};
}  // namespace OpenABF