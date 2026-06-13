/*
OpenABF
https://gitlab.com/educelab/OpenABF

Copyright 2025 EduceLab

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/
#pragma once

// clang-format off
// #include "OpenABF/Exceptions.hpp"


#include <stdexcept>
#include <string>

namespace OpenABF
{

/** @brief Solver exception */
class SolverException : public std::runtime_error
{
public:
    /** @brief Constructor with message */
    explicit SolverException(const char* msg) : std::runtime_error(msg) {}
    /** @brief Constructor with message */
    explicit SolverException(const std::string& msg) : std::runtime_error(msg) {}
};

/** @brief Solver exception */
class MeshException : public std::runtime_error
{
public:
    /** @brief Constructor with message */
    explicit MeshException(const char* msg) : std::runtime_error(msg) {}
    /** @brief Constructor with message */
    explicit MeshException(const std::string& msg) : std::runtime_error(msg) {}
};

}  // namespace OpenABF
// #include "OpenABF/Math.hpp"


#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace OpenABF
{

/** @brief Pi, templated for floating-point type */
template <class T>
constexpr T PI = T(3.1415926535897932385L);

/** @brief Inf, templated for floating-point type */
template <class T>
constexpr T INF = std::numeric_limits<T>::infinity();

/** @brief Vector dot product (inner product) */
template <typename T1, typename T2>
auto dot(const T1& a, const T2& b)
{
    using Ret = decltype(*std::begin(a));
    return std::inner_product(std::begin(a), std::end(a), std::begin(b), Ret(0));
}
/** @brief Vector cross product */
template <typename T1, typename T2>
auto cross(const T1& a, const T2& b) -> T1
{
    T1 c;
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
    return c;
}

/** @brief Norm type enumeration */
enum class Norm {
    L1,  /**< \f$L_1\f$ norm */
    L2,  /**< \f$L_2\f$ norm */
    LInf /**< \f$L_{Inf}\f$ norm */
};

/** @brief Compute vector norm */
template <class Vector>
auto norm(const Vector& v, Norm norm = Norm::L2)
{
    using Ret = decltype(*std::begin(v));
    switch (norm) {
        case Norm::L1: {
            return std::accumulate(std::begin(v), std::end(v), Ret(0),
                                   [](auto a, auto b) { return a + std::abs(b); });
        }
        case Norm::L2: {
            auto sum = std::accumulate(std::begin(v), std::end(v), Ret(0),
                                       [](auto a, auto b) { return a + (b * b); });
            return std::sqrt(sum);
        }
        case Norm::LInf: {
            return std::abs(*std::max_element(std::begin(v), std::end(v), [](auto a, auto b) {
                return std::abs(a) < std::abs(b);
            }));
        }
    }
    throw std::invalid_argument("Invalid norm option");
}

/** @brief Normalize a vector (i.e. compute a unit vector) */
template <class Vector>
auto normalize(Vector v)
{
    return v / norm(v, Norm::L2);
}

/** @brief Compute the interior angle between two vectors */
template <class Vector1, class Vector2>
auto interior_angle(const Vector1& a, const Vector2& b)
{
    return std::acos(dot(a, b) / (norm(a, Norm::L2) * norm(b, Norm::L2)));
}

/** @brief Convert degrees to radians */
template <typename T = float, typename T2,
          std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
constexpr auto to_radians(T2 deg) -> T
{
    return deg * PI<T> / T(180);
}

/** @brief Convert radians to degrees */
template <typename T = float, typename T2,
          std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
constexpr auto to_degrees(T2 rad) -> T
{
    return rad * T(180) / PI<T>;
}

}  // namespace OpenABF
// #include "OpenABF/Vec.hpp"


#include <array>
#include <iostream>

// #include "OpenABF/Math.hpp"


namespace OpenABF
{
/**
 * @brief N-dimensional vector class
 *
 * Essentially a wrapper around std::array that makes it more convenient for
 * vector math purposes.
 *
 * @tparam T Element type
 * @tparam Dims Number of elements
 */
template <typename T, std::size_t Dims, std::enable_if_t<std::is_arithmetic<T>::value, bool> = true>
class Vec
{
    /** Underlying element storage */
    using Container = std::array<T, Dims>;

public:
    /** Element type */
    using value_type = T;
    /** Vector size type */
    using size_type = typename Container::size_type;
    /** Difference type */
    using difference_type = typename Container::difference_type;
    /** Reference type */
    using reference = value_type&;
    /** Const reference type */
    using const_reference = const value_type&;
    /** Pointer type */
    using pointer = value_type*;
    /** Const pointer type */
    using const_pointer = const value_type*;
    /** Iterator type */
    using iterator = typename Container::iterator;
    /** Const iterator type */
    using const_iterator = typename Container::const_iterator;
    /** Reverse iterator type */
    using reverse_iterator = typename Container::reverse_iterator;
    /** Const reverse iterator type */
    using const_reverse_iterator = typename Container::const_reverse_iterator;

    /** @brief Default constructor */
    Vec() { val_.fill(0); }

    /**
     * @brief Construct with element values
     *
     * The number of arguments provided must match Dims.
     */
    template <typename... Args>
    explicit Vec(Args... args)
    {
        static_assert(sizeof...(args) == Dims, "Incorrect number of arguments");
        std::size_t i{0};
        ((val_[i++] = args), ...);
    }

    /** @brief Copy constructor */
    template <typename Vector>
    explicit Vec(const Vector& vec)
    {
        std::copy(std::begin(vec), std::end(vec), val_.begin());
    }

    /** @brief Bounds-checked element access */
    constexpr reference at(size_type pos) { return val_.at(pos); }
    /** @brief Bounds-checked element access */
    constexpr const_reference at(size_type pos) const { return val_.at(pos); }
    /** @brief Element access */
    constexpr reference operator[](size_type i) { return val_[i]; }
    /** @brief Element access */
    constexpr const_reference operator[](size_type i) const { return val_[i]; }

    /** @brief First element */
    constexpr reference front() { return val_.front(); }
    /** @brief First element */
    constexpr const_reference front() const { return val_.front(); }
    /** @brief Last element */
    constexpr reference back() { return val_.back(); }
    /** @brief Last element */
    constexpr const_reference back() const { return val_.back(); }

    /** @brief Get a pointer to the first element of the raw data */
    constexpr pointer data() { return val_.data(); }
    /** @brief Get a pointer to the first element of the raw data */
    constexpr const_pointer data() const { return val_.data(); }

    /** @brief Get an iterator to the first element of the vector */
    constexpr iterator begin() noexcept { return val_.begin(); }
    /** @brief Get an iterator to the first element of the vector */
    constexpr const_iterator begin() const noexcept { return val_.begin(); }
    /** @brief Get an iterator to the first element of the vector */
    constexpr const_iterator cbegin() const noexcept { return val_.cbegin(); }

    /** @brief Get an iterator to one past the last element in the vector */
    constexpr iterator end() noexcept { return val_.end(); }
    /** @brief Get an iterator to one past the last element in the vector */
    constexpr const_iterator end() const noexcept { return val_.end(); }
    /** @brief Get an iterator to one past the last element in the vector */
    constexpr const_iterator cend() const noexcept { return val_.cend(); }

    /** @brief Get an iterator to the first element of the reverse vector */
    constexpr reverse_iterator rbegin() noexcept { return val_.rbegin(); }
    /** @brief Get an iterator to the first element of the reverse vector */
    constexpr const_reverse_iterator rbegin() const noexcept { return val_.rbegin(); }
    /** @brief Get an iterator to the first element of the reverse vector */
    constexpr const_reverse_iterator crbegin() const noexcept { return val_.crbegin(); }

    /**
     * @brief Get an iterator to one past the last element in the reverse vector
     */
    constexpr reverse_iterator rend() noexcept { return val_.rend(); }
    /**
     * @brief Get an iterator to one past the last element in the reverse vector
     */
    constexpr const_reverse_iterator rend() const noexcept { return val_.rend(); }
    /**
     * @brief Get an iterator to one past the last element in the reverse vector
     */
    constexpr const_reverse_iterator crend() const noexcept { return val_.crend(); }

    /** @brief Return whether the vector is empty (uninitialized) */
    constexpr bool empty() const noexcept { return val_.empty(); }
    /** @brief Return the number of elements in the vector */
    constexpr size_type size() const noexcept { return val_.size(); }

    /** @brief Fill the vector with a value */
    constexpr void fill(const T& value) { val_.fill(value); }
    /** @brief Swap this vector with another vector */
    constexpr void swap(Vec& other) noexcept { val_.swap(other.val_); }

    /** @brief Equality comparison operator */
    bool operator==(const Vec& rhs) const { return val_ == rhs.val_; }
    /** @brief Inequality comparison operator */
    bool operator!=(const Vec& rhs) const { return val_ != rhs.val_; }

    /** @brief Assignment operator */
    template <class Vector>
    Vec& operator=(const Vector& b)
    {
        std::size_t idx{0};
        for (auto& v : val_) {
            v = b[idx++];
        }
        return *this;
    }

    /** @brief Assignment operator for std::initializer_list */
    template <typename T2>
    Vec& operator=(const std::initializer_list<T2>& b)
    {
        auto it = b.begin();
        for (auto& v : val_) {
            v = *it;
            it++;
        }
        return *this;
    }

    /** @brief Addition assignment operator */
    template <class Vector>
    Vec& operator+=(const Vector& b)
    {
        std::size_t idx{0};
        for (auto& v : val_) {
            v += b[idx++];
        }
        return *this;
    }

    /** @brief Addition operator */
    template <class Vector>
    friend Vec operator+(Vec lhs, const Vector& rhs)
    {
        lhs += rhs;
        return lhs;
    }

    /** @brief Subtraction assignment operator */
    template <class Vector>
    Vec& operator-=(const Vector& b)
    {
        std::size_t idx{0};
        for (auto& v : val_) {
            v -= b[idx++];
        }
        return *this;
    }

    /** @brief Subtraction operator */
    template <class Vector>
    friend Vec operator-(Vec lhs, const Vector& rhs)
    {
        lhs -= rhs;
        return lhs;
    }

    /** @brief Multiplication assignment operator */
    template <typename T2, std::enable_if_t<std::is_arithmetic<T2>::value, bool> = true>
    Vec& operator*=(const T2& b)
    {
        for (auto& v : val_) {
            v *= b;
        }
        return *this;
    }

    /** @brief Multiplication operator */
    template <typename T2, std::enable_if_t<std::is_arithmetic<T2>::value, bool> = true>
    friend Vec operator*(Vec lhs, const T2& rhs)
    {
        lhs *= rhs;
        return lhs;
    }

    /** @brief Division assignment operator */
    template <typename T2, std::enable_if_t<std::is_arithmetic<T2>::value, bool> = true>
    Vec& operator/=(const T2& b)
    {
        for (auto& v : val_) {
            v /= b;
        }
        return *this;
    }

    /** @brief Division operator */
    template <typename T2, std::enable_if_t<std::is_arithmetic<T2>::value, bool> = true>
    friend Vec operator/(Vec lhs, const T2& rhs)
    {
        lhs /= rhs;
        return lhs;
    }

    /** @brief Compute the vector dot product (i.e. inner product) */
    template <class Vector>
    T dot(const Vector& v)
    {
        return OpenABF::dot(val_, v);
    }

    /** @brief Compute the vector cross product */
    template <class Vector, std::size_t D = Dims>
    std::enable_if_t<D == 3, Vec> cross(const Vector& v)
    {
        return OpenABF::cross(*this, v);
    }

    /** @brief Compute the vector magnitude */
    T magnitude() const { return OpenABF::norm(*this, Norm::L2); }

    /** @brief Return the unit vector of this vector */
    Vec unit() const { return OpenABF::normalize(*this); }

private:
    /** Values */
    Container val_{};
};

/** @brief 3D, 32-bit float vector */
using Vec3f = Vec<float, 3>;
/** @brief 3D, 64-bit float vector */
using Vec3d = Vec<double, 3>;

/** Debug: Print a vector to a std::ostream */
template <typename T, std::size_t Dims>
std::ostream& operator<<(std::ostream& os, const Vec<T, Dims>& vec)
{
    os << "[";
    std::size_t i{0};
    for (const auto& v : vec) {
        if (i++ > 0) {
            os << ", ";
        }
        os << v;
    }
    os << "]";
    return os;
}

}  // namespace OpenABF


// #include "OpenABF/HalfEdgeMesh.hpp"


#include <algorithm>
#include <array>
#include <iterator>
#include <map>
#include <memory>
#include <optional>
#include <queue>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// #include "OpenABF/Exceptions.hpp"

// #include "OpenABF/Vec.hpp"


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
    /** First iterator in the range. */
    Iter first_;
    /** One-past-last iterator in the range. */
    Iter last_;
    /** Range begin iterator. */
    auto begin() const -> Iter { return first_; }
    /** Range end iterator. */
    auto end() const -> Iter { return last_; }
    /** True if the range is empty. */
    auto empty() const -> bool { return first_ == last_; }
    /** First element of the range. */
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
    /** Difference type. */
    using difference_type = std::ptrdiff_t;
    /** Value type. */
    using value_type = typename std::iterator_traits<Iter>::value_type;
    /** Pointer type. */
    using pointer = typename std::iterator_traits<Iter>::pointer;
    /** Reference type. */
    using reference = typename std::iterator_traits<Iter>::reference;
    /** Iterator category. */
    using iterator_category = std::input_iterator_tag;

    /** Default constructor. */
    FilteringIterator() = default;
    /** Construct over `[current, end)` filtering by `pred`. */
    FilteringIterator(Iter current, Iter end, Pred pred) : current_{current}, end_{end}, pred_{pred}
    {
        advance_to_next();
    }

    /** Dereference. */
    auto operator*() const -> reference { return *current_; }
    /** Member access. */
    auto operator->() const -> pointer { return &*current_; }

    /** Pre-increment; skips elements failing `pred`. */
    auto operator++() -> FilteringIterator&
    {
        ++current_;
        advance_to_next();
        return *this;
    }

    /** Equality. */
    auto operator==(const FilteringIterator& other) const -> bool
    {
        return current_ == other.current_;
    }
    /** Inequality. */
    auto operator!=(const FilteringIterator& other) const -> bool { return !(*this == other); }

private:
    /** Advance `current_` past elements failing `pred_`. */
    void advance_to_next()
    {
        while (current_ != end_ && !pred_(*current_)) {
            ++current_;
        }
    }

    /** Underlying iterator position. */
    Iter current_{};
    /** Underlying end iterator. */
    Iter end_{};
    /** Predicate selecting which elements to visit. */
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

        /** Dereference (const overload) */
        template <bool C = Const>
        auto operator*() const -> std::enable_if_t<C, reference>
        {
            return current_;
        }
        /** Dereference (non-const overload) */
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
        /** Advance `current_` past boundary edges, wrapping to `nullptr` at end. */
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

        /** First edge in the wheel (used to detect wrap-around). */
        EdgePtr head_{};
        /** Current edge in the wheel; `nullptr` at end. */
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

        /** Underlying face-vector iterator type */
        using FaceVecIter = typename std::vector<FacePtr>::const_iterator;

        /** Construct at position (begin or end depending on faceIt == faceEnd) */
        EdgesIterator(FaceVecIter faceIt, FaceVecIter faceEnd) : faceIt_{faceIt}, faceEnd_{faceEnd}
        {
            if (faceIt_ != faceEnd_) {
                edgeIt_ = FaceIterator<true>{(*faceIt_)->head, (*faceIt_)->head};
                advance_if_face_exhausted();
            }
        }

        /** Dereference (const overload) */
        template <bool C = Const>
        auto operator*() const -> std::enable_if_t<C, reference>
        {
            return *edgeIt_;
        }
        /** Dereference (non-const overload) */
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
        /** Advance `faceIt_` until a face with remaining edges is found, or the end is reached. */
        void advance_if_face_exhausted()
        {
            while (edgeIt_ == FaceIterator<true>() && faceIt_ != faceEnd_) {
                ++faceIt_;
                if (faceIt_ != faceEnd_) {
                    edgeIt_ = FaceIterator<true>{(*faceIt_)->head, (*faceIt_)->head};
                }
            }
        }

        /** Current position in the face vector. */
        FaceVecIter faceIt_{};
        /** End of the face vector. */
        FaceVecIter faceEnd_{};
        /** Edge iterator within the current face. */
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
     * @brief Insert new vertices from a brace-enclosed initializer list
     *
     * Convenience overload of insert_vertices() accepting nested
     * `std::initializer_list` syntax, e.g. `{{0,0,0}, {1,0,0}, {0,1,0}}`.
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

    /**
     * @brief One element of `extract_connected_components()`'s output.
     *
     * Holds an extracted sub-mesh and the two index maps callers need to
     * relate sub-mesh elements back to the source mesh. The mapping is 1-to-1
     * in both directions for a manifold input — every source vertex/face is
     * in exactly one component.
     */
    struct ExtractedComponent {
        /** Deep copy of the connected component as an independent mesh. */
        Pointer mesh;
        /** `vertex_map[sub_idx] == original_idx`. */
        std::vector<std::size_t> vertex_map;
        /** `face_map[sub_idx] == original_idx`. */
        std::vector<std::size_t> face_map;
    };

    /**
     * @brief Extract each connected component of this mesh into its own
     * independent mesh.
     *
     * Returns one `ExtractedComponent` per connected component. The extracted
     * meshes are deep copies that preserve VertexTraits, EdgeTraits, and
     * FaceTraits via the relevant copy constructors. Vertex and face indices
     * in the extracted meshes are re-densified to `0..N-1`; the `vertex_map`
     * and `face_map` let callers scatter/gather per-vertex and per-face data
     * between the sub-mesh and the source.
     *
     * A single-component mesh produces a vector of length 1 whose extracted
     * mesh is equivalent to `clone()` and whose maps are both the identity.
     */
    auto extract_connected_components() -> std::vector<ExtractedComponent>
    {
        std::vector<ExtractedComponent> result;
        for (const auto& component : connected_components()) {
            ExtractedComponent out;
            out.mesh = HalfEdgeMesh::New();
            std::unordered_map<std::size_t, std::size_t> remap;
            out.face_map.reserve(component.size());

            // Insert vertices first so the sub-mesh has the targets that
            // clone_face_ will look up via verts_.at(...). Copy ctor here
            // preserves VertexTraits.
            for (const auto& face : component) {
                for (const auto& edge : *face) {
                    const auto& v = edge->vertex;
                    if (remap.find(v->idx) == remap.end()) {
                        const auto newIdx = out.mesh->insert_vertex(*v);
                        out.mesh->verts_[newIdx]->edge = nullptr;
                        remap.emplace(v->idx, newIdx);
                        out.vertex_map.push_back(v->idx);
                    }
                }
            }

            // Clone each face into the sub-mesh, remapping source vertex
            // indices to the sub-mesh's densified indices. clone_face_
            // preserves FaceTraits and EdgeTraits. Face indices in the sub-
            // mesh come out densified in insertion order, so the loop's idx
            // matches the sub-face idx and we can record face_map alongside.
            const auto remapFn = [&remap](std::size_t i) { return remap.at(i); };
            for (const auto& face : component) {
                out.mesh->clone_face_(face, remapFn);
                out.face_map.push_back(face->idx);
            }
            out.mesh->update_boundary();

            result.push_back(std::move(out));
        }
        return result;
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
        rekey_edge_to_vertex_(newFwd->next, newEnd);

        // Update new boundary edges' next/prev. The rekey_* calls move
        // half-edges between edges_ multimap buckets so that
        // outgoing_edges(idx) / Vertex::is_manifold() remain consistent.
        if (startOnBoundary) {
            rekey_edge_to_vertex_(startOut, newStart);
            newBwd->next = startOut;
            startOut->prev = newBwd;
            oldFwd->prev = startIn;
            startIn->next = oldFwd;
            for (auto e : newStart->wheel()) {
                rekey_edge_to_vertex_(e, newStart);
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
                rekey_edge_to_vertex_(e, newEnd);
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
     * Move @p e from its current edges_ bucket to the bucket keyed by
     * @p newVert, and update e->vertex.
     *
     * edges_ is a multimap keyed by half-edge origin vertex index. Any code
     * that reassigns a half-edge's origin (e.g. split_edge duplicating a
     * vertex and moving a fan of half-edges to the new vertex) must call
     * this to keep the multimap consistent, otherwise outgoing_edges() /
     * incoming_edges() / Vertex::is_manifold() return stale results.
     */
    void rekey_edge_to_vertex_(const EdgePtr& e, const VertPtr& newVert)
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
     * from an existing mesh. Preserves FaceTraits (via `Face::New(*face)`) and
     * EdgeTraits (via `Edge::New(*e)`) by exploiting the trait-only copy
     * constructors.
     *
     * @param face Existing face from the mesh being cloned (may belong to a
     *     different mesh than `this`).
     * @param remap Callable mapping source-mesh vertex idx -> this-mesh
     *     vertex idx. The default identity remap is used by `clone()` where
     *     the two meshes share vertex indexing; `extract_connected_components`
     *     passes a real remap when sub-meshes have re-densified indices.
     */
    template <typename RemapFn>
    auto clone_face_(const FacePtr& face, RemapFn remap)
    {
        // Copy the existing face
        auto f = Face::New(*face);

        // Pre-make all edges
        std::vector<std::size_t> idxs;
        for (const auto& e : *face) {
            auto startIdx = remap(e->vertex->idx);
            auto endIdx = remap(e->pair->vertex->idx);
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

    /** Identity-remap overload (source and target share vertex indexing). */
    auto clone_face_(const FacePtr& face)
    {
        return clone_face_(face, [](std::size_t i) { return i; });
    }
};
}  // namespace OpenABF

// #include "OpenABF/ABF.hpp"


#include <cassert>
#include <cmath>
#include <limits>
#include <vector>

#include <Eigen/SparseLU>

// #include "OpenABF/Exceptions.hpp"

// #include "OpenABF/HalfEdgeMesh.hpp"

// #include "OpenABF/Math.hpp"


namespace OpenABF
{

namespace traits
{
/** @brief ABF and ABFPlusPlus vertex traits */
template <typename T>
struct ABFVertexTraits : DefaultVertexTraits<T> {
    /** Lagrange Multiplier: Planarity constraint */
    T lambda_plan{0};
    /** Lagrange Multiplier: Reconstruction constraint */
    T lambda_len{1};
};

/** @brief ABF and ABFPlusPlus edge traits */
template <typename T>
struct ABFEdgeTraits : DefaultEdgeTraits<T> {
    /** 3D angle */
    T beta{0};
    /** Optimal (i.e. target) angle */
    T phi{0};
    /** Angle weight */
    T weight{0};
    /** Sin of alpha, because it's used a lot */
    T alpha_sin{0};
    /** Cos of alpha, because it's used a lot */
    T alpha_cos{0};
};

/** @brief ABF and ABFPlusPlus face traits */
template <typename T>
struct ABFFaceTraits : DefaultFaceTraits<T> {
    /** Lagrange Multiplier: Triangle validity constraint */
    T lambda_tri{0};
};
}  // namespace traits

/** @brief %ABF and ABF++ implementation details */
namespace detail::ABF
{

/** @brief A HalfEdgeMesh with the %ABF traits */
template <typename T>
using Mesh = HalfEdgeMesh<T, 3, traits::ABFVertexTraits<T>, traits::ABFEdgeTraits<T>,
                          traits::ABFFaceTraits<T>>;

/** @brief Initialize the %ABF angles and weights from the edge alpha values */
template <typename T, class MeshPtr>
void InitializeAnglesAndWeights(MeshPtr& m)
{
    // Initialize and bound angle properties
    static constexpr auto MinAngle = PI<T> / T(180);
    static constexpr auto MaxAngle = PI<T> - MinAngle;
    for (auto& e : m->edges()) {
        e->alpha = e->beta = e->phi = std::min(std::max(e->alpha, MinAngle), MaxAngle);
        e->alpha_sin = std::sin(e->alpha);
        e->alpha_cos = std::cos(e->alpha);
        e->weight = T(1) / (e->phi * e->phi);
    }

    // Update weights for interior vertices
    for (auto& v : m->vertices_interior()) {
        auto wheel = v->wheel();
        auto angle_sum = std::accumulate(wheel.begin(), wheel.end(), T(0),
                                         [](auto a, auto b) { return a + b->beta; });
        for (auto& e : wheel) {
            e->phi *= 2 * PI<T> / angle_sum;
            e->weight = T(1) / (e->phi * e->phi);
        }
    }
}

/** @brief Compute ∇CTri w.r.t LambdaTri == CTri */
template <typename T, class FacePtr, std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
auto TriGrad(const FacePtr& f) -> T
{
    T g = -PI<T>;
    for (const auto& e : *f) {
        g += e->alpha;
    }
    return g;
}

/** @brief Compute ∇CPlan w.r.t LambdaPlan == CPlan */
template <typename T, class VertPtr, std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
auto PlanGrad(const VertPtr& v) -> T
{
    T g = -2 * PI<T>;
    for (const auto& e : v->wheel()) {
        g += e->alpha;
    }
    return g;
}

/** @brief Compute ∇CLen w.r.t LambdaLen == CLen */
template <typename T, class VertPtr, std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
auto LenGrad(const VertPtr& vertex) -> T
{
    T p1{1};
    T p2{1};
    for (const auto& e : vertex->wheel()) {
        p1 *= e->next->alpha_sin;
        p2 *= e->next->next->alpha_sin;
    }
    return p1 - p2;
}

/** @brief Compute ∇CLen w.r.t edge->alpha */
template <typename T, class VertPtr, class EdgePtr,
          std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
auto LenGrad(const VertPtr& vertex, const EdgePtr& edge) -> T
{
    T p1{1};
    T p2{1};
    for (const auto& a : vertex->wheel()) {
        auto b = a->next;
        if (b == edge) {
            p1 *= b->alpha_cos;
            p2 = T(0);
        } else {
            p1 *= b->alpha_sin;
        }

        auto c = a->next->next;
        if (c == edge) {
            p1 = T(0);
            p2 *= c->alpha_cos;
        } else {
            p2 *= c->alpha_sin;
        }
    }
    return p1 - p2;
}

/** @brief Compute ∇F w.r.t an edge's alpha */
template <typename T, class EdgePtr, std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
auto AlphaGrad(const EdgePtr& edge) -> T
{
    // δE/δα
    auto g = (edge->alpha - edge->phi) * edge->weight;
    // δCTri/δα
    g += edge->face->lambda_tri;
    for (const auto& e : *edge->face) {
        // Skip boundary vertices
        if (e->vertex->is_boundary()) {
            continue;
        }
        if (e == edge) {
            // δCPlan/δα
            g += e->vertex->lambda_plan;
        } else {
            // δCLen/δα
            auto d = LenGrad<T>(e->vertex, edge);
            d *= e->vertex->lambda_len;
            g += d;
        }
    }
    return g;
}

/** @brief Compute ∇F w.r.t all parameters */
template <typename T, class MeshPtr, std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
auto Gradient(const MeshPtr& mesh) -> T
{
    T g{0};
    for (const auto& f : mesh->faces()) {
        // AlphaGrad for all edges
        for (const auto& e : *f) {
            auto gAlpha = AlphaGrad<T>(e);
            g += gAlpha * gAlpha;
        }
        // TriGrad for all faces
        auto gTri = TriGrad<T>(f);
        g += gTri * gTri;
    }

    // PlanGrad and LenGrad for all interior vertices
    for (const auto& v : mesh->vertices_interior()) {
        auto gPlan = PlanGrad<T>(v);
        g += gPlan * gPlan;

        auto gLen = LenGrad<T>(v);
        g += gLen * gLen;
    }
    return g;
}
}  // namespace detail::ABF

/**
 * @brief Compute parameterized interior angles using Angle-based flattening
 *
 * Iteratively computes a new set of interior angles which minimize the total
 * angular error of the parameterized mesh. This follows the standard
 * angled-based flattening formulation, which directly solves for the objective
 * functions and constraints. ABFPlusPlus is generally preferred as it
 * dramatically simplifies the size of the solved problem without introducing
 * more error.
 *
 * This class **does not** compute a parameterized mesh. Rather, it calculates
 * the optimal interior angles for such a mesh. To convert this information
 * into a full parameterization, pass the processed HalfEdgeMesh to
 * AngleBasedLSCM.
 *
 * Implements "Parameterization of faceted surfaces for meshing using
 * angle-based flattening" by Sheffer and de Sturler (2001)
 * \cite sheffer2001abf.
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
class ABF
{
public:
    /** @brief Mesh type alias */
    using Mesh = MeshType;

    /** @brief Set the maximum number of iterations */
    void setMaxIterations(const std::size_t it) { maxIters_ = it; }

    /** @brief Set the gradient convergence threshold */
    void setGradientThreshold(T t) { gradThreshold_ = t; }

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

    /** @copydoc ABF::Compute */
    void compute(typename Mesh::Pointer& mesh)
    {
        Compute(mesh, iters_, grad_, maxIters_, gradThreshold_);
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
        auto vCnt = mesh->num_vertices();
        std::vector<std::size_t> vIdx2vIntIdx(vCnt, std::numeric_limits<std::size_t>::max());
        {
            std::size_t newIdx{0};
            for (const auto& v : mesh->vertices_interior()) {
                vIdx2vIntIdx[v->idx] = newIdx++;
            }
        }

        while (gradient > gradThreshold and gradDelta > gradThreshold and iters < maxIters) {
            if (std::isnan(gradient) or std::isinf(gradient)) {
                throw MeshException("Mesh gradient cannot be computed");
            }
            // Typedefs
            using Triplet = Eigen::Triplet<T>;
            using SparseMatrix = Eigen::SparseMatrix<T>;
            using DenseVector = Eigen::Matrix<T, Eigen::Dynamic, 1>;

            // Helpful parameters
            auto vIntCnt = mesh->num_vertices_interior();
            auto edgeCnt = mesh->num_edges();
            auto faceCnt = mesh->num_faces();

            //// RHS ////
            // b1 = -alpha gradient
            std::vector<Triplet> triplets;
            std::size_t idx{0};
            for (const auto& e : mesh->edges()) {
                triplets.emplace_back(idx, 0, -AlphaGrad<T>(e));
                ++idx;
            }

            // b2 = -lambda gradient
            // lambda tri
            for (const auto& f : mesh->faces()) {
                triplets.emplace_back(idx, 0, -TriGrad<T>(f));
                ++idx;
            }
            // lambda plan and lambda len
            for (const auto& v : mesh->vertices_interior()) {
                triplets.emplace_back(idx, 0, -PlanGrad<T>(v));
                triplets.emplace_back(vIntCnt + idx, 0, -LenGrad<T>(v));
                ++idx;
            }
            SparseMatrix b(edgeCnt + faceCnt + 2 * vIntCnt, 1);
            b.reserve(triplets.size());
            b.setFromTriplets(triplets.begin(), triplets.end());

            ///// LHS /////
            // Lambda = diag(2/w)
            // v.weight == 1/w, so Lambda is diag(2*weight)
            // We only need Lambda Inverse, so this is 1 / 2*weight
            triplets.clear();
            idx = 0;
            for (const auto& e : mesh->edges()) {
                triplets.emplace_back(idx, idx, 2 * e->weight);
                ++idx;
            }

            // J
            // Jacobian of the CTri constraints
            for (idx = 0; idx < faceCnt; idx++) {
                auto row = idx + edgeCnt;
                auto col = 3 * idx;
                triplets.emplace_back(row, col, 1);
                triplets.emplace_back(row, col + 1, 1);
                triplets.emplace_back(row, col + 2, 1);

                // Jt
                triplets.emplace_back(col, row, 1);
                triplets.emplace_back(col + 1, row, 1);
                triplets.emplace_back(col + 2, row, 1);
            }
            for (const auto& v : mesh->vertices_interior()) {
                auto row = idx + edgeCnt;
                for (const auto& e0 : v->wheel()) {
                    // Jacobian of the CPlan constraint
                    triplets.emplace_back(row, e0->idxI.value(), 1);
                    triplets.emplace_back(e0->idxI.value(), row, 1);

                    // Jacobian of the CLen constraint
                    auto e1 = e0->next;
                    auto e2 = e1->next;
                    auto d1 = LenGrad<T>(v, e1);
                    auto d2 = LenGrad<T>(v, e2);
                    triplets.emplace_back(vIntCnt + row, e1->idxI.value(), d1);
                    triplets.emplace_back(vIntCnt + row, e2->idxI.value(), d2);
                    triplets.emplace_back(e1->idxI.value(), vIntCnt + row, d1);
                    triplets.emplace_back(e2->idxI.value(), vIntCnt + row, d2);
                }
                ++idx;
            }
            auto Asize = edgeCnt + faceCnt + 2 * vIntCnt;
            SparseMatrix A(Asize, Asize);
            A.reserve(triplets.size());
            A.setFromTriplets(triplets.begin(), triplets.end());

            A.makeCompressed();
            Solver solver;
            solver.compute(A);
            if (solver.info() != Eigen::ComputationInfo::Success) {
                throw SolverException("ABF: Failed to solve A");
            }
            DenseVector delta = solver.solve(b);
            if (solver.info() != Eigen::ComputationInfo::Success) {
                throw SolverException("ABF: Failed to solve b");
            }

            // alpha += delta_alpha
            // Update sin and cos
            idx = 0;
            for (auto& e : mesh->edges()) {
                e->alpha += delta(idx++, 0);
                e->alpha = std::min(std::max(e->alpha, T(0)), PI<T>);
                e->alpha_sin = std::sin(e->alpha);
                e->alpha_cos = std::cos(e->alpha);
            }

            // lambda += delta_lambda
            for (auto& f : mesh->faces()) {
                f->lambda_tri += delta(idx++, 0);
            }
            auto base = edgeCnt + faceCnt;
            for (auto& v : mesh->vertices_interior()) {
                auto intIdx = vIdx2vIntIdx[v->idx];
                assert(intIdx != std::numeric_limits<std::size_t>::max());
                v->lambda_plan += delta(base + intIdx, 0);
                v->lambda_len += delta(base + vIntCnt + intIdx, 0);
            }

            // Recalculate gradient for next iteration
            auto newGrad = detail::ABF::Gradient<T>(mesh);
            gradDelta = std::abs(newGrad - gradient);
            gradient = newGrad;
            iters++;
        }
    }

    /** @copydoc ABF::Compute */
    static void Compute(typename Mesh::Pointer& mesh)
    {
        std::size_t iters{0};
        T grad{0};
        Compute(mesh, iters, grad);
    }

protected:
    /** Gradient */
    T grad_{0};
    /** Number of executed iterations */
    std::size_t iters_{0};
    /** Max iterations */
    std::size_t maxIters_{10};
    /** Gradient convergence threshold */
    T gradThreshold_{0.001};
};

}  // namespace OpenABF
// #include "OpenABF/ABFPlusPlus.hpp"


#include <cassert>
#include <cmath>
#include <limits>
#include <vector>

#include <Eigen/SparseLU>

// #include "OpenABF/ABF.hpp"

// #include "OpenABF/Exceptions.hpp"

// #include "OpenABF/HalfEdgeMesh.hpp"

// #include "OpenABF/Math.hpp"


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
    void setMaxIterations(std::size_t it) { maxIters_ = it; }

    /** @brief Set the gradient convergence threshold */
    void setGradientThreshold(T t) { gradThreshold_ = t; }

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
        Compute(mesh, iters_, grad_, maxIters_, gradThreshold_);
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

        while (gradient > gradThreshold and gradDelta > gradThreshold and iters < maxIters) {
            if (std::isnan(gradient) or std::isinf(gradient)) {
                throw MeshException("Mesh gradient cannot be computed");
            }
            // Typedefs
            using Triplet = Eigen::Triplet<T>;
            using SparseMatrix = Eigen::SparseMatrix<T>;
            using DenseVector = Eigen::Matrix<T, Eigen::Dynamic, 1>;

            // Helpful parameters
            auto vIntCnt = mesh->num_vertices_interior();
            auto edgeCnt = mesh->num_edges();
            auto faceCnt = mesh->num_faces();

            // b1 = -alpha gradient
            std::vector<Triplet> triplets;
            std::size_t idx{0};
            for (const auto& e : mesh->edges()) {
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
            SparseMatrix b2(faceCnt + 2 * vIntCnt, 1);
            b2.reserve(triplets.size());
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
            SparseMatrix J(faceCnt + 2 * vIntCnt, 3 * faceCnt);
            J.reserve(triplets.size());
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
            SparseMatrix LambdaInv(edgeCnt, edgeCnt);
            LambdaInv.reserve(edgeCnt);
            LambdaInv.setFromTriplets(triplets.begin(), triplets.end());

            // solve Eq. 16
            auto bstar = J * LambdaInv * b1 - b2;
            auto JLiJt = J * LambdaInv * J.transpose();

            SparseMatrix LambdaStarInv = JLiJt.block(0, 0, faceCnt, faceCnt);
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
            SparseMatrix A = Jstar * LambdaStarInv * JstarT - Jstar2;
            SparseMatrix b = Jstar * LambdaStarInv * bstar1 - bstar2;
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
            DenseVector deltaLambda(deltaLambda1.rows() + deltaLambda2.rows(), 1);
            deltaLambda << DenseVector(deltaLambda1), DenseVector(deltaLambda2);

            // Compute Eq. 10 -> delta_alpha
            DenseVector deltaAlpha = LambdaInv * (b1 - J.transpose() * deltaLambda);

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
    std::size_t maxIters_{10};
    /** Gradient convergence threshold */
    T gradThreshold_{0.001};
};

}  // namespace OpenABF
// #include "OpenABF/AngleBasedLSCM.hpp"


#include <optional>
#include <type_traits>
#include <utility>

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>

// #include "OpenABF/Exceptions.hpp"

// #include "OpenABF/HalfEdgeMesh.hpp"

// #include "OpenABF/Math.hpp"

// #include "OpenABF/detail/LSCMSystem.hpp"


#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <unordered_map>
#include <vector>

#include <Eigen/SparseCore>

namespace OpenABF::detail::lscm
{

/**
 * @brief Outputs of `buildSystem`: the LSCM least-squares system for a mesh
 *        with two pinned vertices.
 *
 * @tparam T Floating-point type
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
 * @brief Build the LSCM sparse system for a mesh with two pinned vertices.
 *
 * Mutates the mesh: places `p0` at the UV origin and `p1` on whichever XY
 * axis its displacement from `p0` has the largest magnitude — same pin
 * placement convention used by `AngleBasedLSCM::ComputeImpl` and
 * `HierarchicalLSCM::solveLSCMLevel`.
 *
 * Assembly follows Lévy et al. 2002 Eq. 10 using the per-edge `alpha` angles
 * already stored on the mesh.
 */
template <typename T, class MeshType>
auto buildSystem(const typename MeshType::Pointer& mesh, const typename MeshType::VertPtr& p0,
                 const typename MeshType::VertPtr& p1) -> SystemParts<T>
{
    using Triplet = Eigen::Triplet<T>;
    using SparseMatrix = Eigen::SparseMatrix<T>;

    // Map selected edge to closest XY axis. Use sign to select direction.
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

    const auto numFaces = mesh->num_faces();
    const auto numVerts = mesh->num_vertices();
    constexpr std::size_t numFixed = 2;
    const auto numFree = numVerts - numFixed;

    // Permutation for free vertices: maps mesh vertex idx → row-pair slot in A.
    std::unordered_map<std::size_t, std::size_t> freeIdxTable;
    freeIdxTable.reserve(numFree);
    for (const auto& v : mesh->vertices()) {
        if (v == p0 or v == p1) {
            continue;
        }
        auto newIdx = freeIdxTable.size();
        freeIdxTable[v->idx] = newIdx;
    }

    // Setup pinned bFixed.
    std::vector<Triplet> tripletsB;
    tripletsB.emplace_back(0, 0, p0->pos[0]);
    tripletsB.emplace_back(1, 0, p0->pos[1]);
    tripletsB.emplace_back(2, 0, p1->pos[0]);
    tripletsB.emplace_back(3, 0, p1->pos[1]);
    SparseMatrix bFixed(2 * numFixed, 1);
    bFixed.reserve(tripletsB.size());
    bFixed.setFromTriplets(tripletsB.begin(), tripletsB.end());

    // Setup variables matrix. Only solving for free vertices, so pins go in
    // a special matrix.
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
// #include "OpenABF/HierarchicalLSCM.hpp"


#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <optional>
#include <queue>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCore>

// #include "OpenABF/AngleBasedLSCM.hpp"

// #include "OpenABF/Exceptions.hpp"

// #include "OpenABF/HalfEdgeMesh.hpp"

// #include "OpenABF/Math.hpp"


namespace OpenABF
{

namespace detail
{
/**
 * @brief Implementation details for HierarchicalLSCM
 */
namespace hlscm
{

/**
 * @brief Symmetric 4×4 quadric matrix for QEM error metric (Garland-Heckbert)
 */
template <typename T>
struct Quadric {
    /** Upper triangle stored row-major: a00 a01 a02 a03 a11 a12 a13 a22 a23 a33 */
    std::array<T, 10> q{};

    Quadric() = default;

    /** Construct from plane equation ax + by + cz + d = 0 */
    Quadric(T a, T b, T c, T d)
        : q{a * a, a * b, a * c, a * d, b * b, b * c, b * d, c * c, c * d, d * d}
    {
    }

    /** In-place accumulation of another quadric */
    auto operator+=(const Quadric& o) -> Quadric&
    {
        for (std::size_t i = 0; i < 10; ++i) {
            q[i] += o.q[i];
        }
        return *this;
    }

    /** Quadric addition */
    friend auto operator+(Quadric a, const Quadric& b) -> Quadric { return a += b; }

    /** Evaluate quadric error at point (x, y, z) */
    auto evaluate(T x, T y, T z) const -> T
    {
        // v^T Q v where Q is the symmetric 4x4 matrix, v = (x, y, z, 1)
        return q[0] * x * x + T(2) * q[1] * x * y + T(2) * q[2] * x * z + T(2) * q[3] * x +
               q[4] * y * y + T(2) * q[5] * y * z + T(2) * q[6] * y + q[7] * z * z +
               T(2) * q[8] * z + q[9];
    }
};

/**
 * @brief Record of a single half-edge collapse for prolongation
 */
template <typename T>
struct CollapseRecord {
    /** Index of the removed vertex (in the original/fine mesh) */
    std::size_t vRemoved;
    /** Index of the kept vertex (in the original/fine mesh) */
    std::size_t vKept;
    /** Post-collapse triangle containing vRemoved (original vertex indices) */
    std::array<std::size_t, 3> containingTri;
    /** Barycentric coordinates of vRemoved in containingTri */
    std::array<T, 3> bary;
};

/**
 * @brief A level in the mesh hierarchy
 */
template <typename T>
struct HierarchyLevel {
    /** Vertex positions (indexed by level-local index) */
    std::vector<Vec<T, 3>> positions;
    /** Face connectivity (each face is 3 level-local indices) */
    std::vector<std::array<std::size_t, 3>> faces;
    /** Map from level-local vertex index to original (finest) vertex index */
    std::vector<std::size_t> localToOriginal;
    /** Map from original vertex index to level-local index. Sized to the
     *  finest-mesh vertex count; vertices not present at this level hold the
     *  sentinel `kAbsent` (`SIZE_MAX`) so callers can distinguish "absent"
     *  from "level-local idx 0".
     */
    std::vector<std::size_t> originalToLocal;

    /** Sentinel marking "no level-local idx for this original vertex". */
    static constexpr std::size_t kAbsent = std::numeric_limits<std::size_t>::max();
};

/**
 * @brief Lightweight flat-array mesh for decimation
 *
 * Copies vertex positions and face connectivity from a HalfEdgeMesh into
 * flat vectors, builds adjacency structures, and supports half-edge collapse.
 */
template <typename T>
class DecimationMesh
{
public:
    /** Build from a HalfEdgeMesh */
    template <class MeshPtr>
    void build(const MeshPtr& mesh, std::size_t pin0, std::size_t pin1)
    {
        auto nv = mesh->num_vertices();
        auto nf = mesh->num_faces();

        positions_.resize(nv);
        alive_.assign(nv, true);
        isBoundary_.assign(nv, false);
        isPinned_.assign(nv, false);
        quadrics_.resize(nv);

        isPinned_[pin0] = true;
        isPinned_[pin1] = true;

        for (const auto& v : mesh->vertices()) {
            positions_[v->idx] = v->pos;
            isBoundary_[v->idx] = v->is_boundary();
        }

        faces_.reserve(nf);
        faceAlive_.reserve(nf);
        vertFaces_.resize(nv);

        for (const auto& f : mesh->faces()) {
            auto e0 = f->head;
            auto e1 = e0->next;
            auto e2 = e1->next;
            std::array<std::size_t, 3> tri{e0->vertex->idx, e1->vertex->idx, e2->vertex->idx};
            auto fi = faces_.size();
            faces_.push_back(tri);
            faceAlive_.push_back(true);
            vertFaces_[tri[0]].push_back(fi);
            vertFaces_[tri[1]].push_back(fi);
            vertFaces_[tri[2]].push_back(fi);
        }

        numAliveVerts_ = nv;
        numAliveFaces_ = nf;

        computeQuadrics_();
        buildEdges_();
    }

    /** Get number of alive vertices */
    [[nodiscard]] auto numAliveVerts() const -> std::size_t { return numAliveVerts_; }

    /** Get number of alive faces */
    [[nodiscard]] auto numAliveFaces() const -> std::size_t { return numAliveFaces_; }

    /**
     * @brief Try to collapse edge (vRemove → vKeep), returning a collapse record
     *
     * Returns nullopt if the collapse is invalid.
     *
     * @param outKeepNbrs  Optional out-vector. On a successful collapse it is
     *                     overwritten with vKeep's post-collapse neighbor
     *                     vertex indices (sorted, unique). Passing a caller-
     *                     owned vector here lets the PQ-update path reuse the
     *                     allocation across collapses instead of letting
     *                     `vertexNeighbors()` allocate a fresh vector each
     *                     time.
     */
    auto tryCollapse(std::size_t vRemove, std::size_t vKeep,
                     std::vector<std::size_t>* outKeepNbrs = nullptr)
        -> std::optional<CollapseRecord<T>>
    {
        if (!alive_[vRemove] || !alive_[vKeep]) {
            return std::nullopt;
        }
        if (isPinned_[vRemove]) {
            return std::nullopt;
        }

        // Find shared faces (will be removed) and vRemove-only faces (will be updated)
        // Reuse scratch storage
        auto& sharedFaces = scratchShared_;
        auto& removeFaces = scratchRemove_;
        sharedFaces.clear();
        removeFaces.clear();

        for (auto fi : vertFaces_[vRemove]) {
            if (!faceAlive_[fi])
                continue;
            bool hasKeep = false;
            for (auto vi : faces_[fi]) {
                if (vi == vKeep) {
                    hasKeep = true;
                    break;
                }
            }
            if (hasKeep)
                sharedFaces.push_back(fi);
            else
                removeFaces.push_back(fi);
        }

        // Interior edges have 2 shared faces; boundary edges have 1.
        if (sharedFaces.empty() || sharedFaces.size() > 2)
            return std::nullopt;

        // Reject collapse of two boundary vertices via an interior edge:
        // vKeep would inherit two disconnected boundary fans → non-manifold.
        if (isBoundary_[vRemove] && isBoundary_[vKeep] && sharedFaces.size() == 2) {
            return std::nullopt;
        }

        // Link condition: collect sorted unique neighbors of vRemove and vKeep
        auto fillSortedNeighbors = [&](std::size_t v, std::vector<std::size_t>& out) {
            out.clear();
            for (auto fi : vertFaces_[v]) {
                if (!faceAlive_[fi])
                    continue;
                for (auto vi : faces_[fi]) {
                    if (vi != v)
                        out.push_back(vi);
                }
            }
            std::sort(out.begin(), out.end());
            out.erase(std::unique(out.begin(), out.end()), out.end());
        };

        fillSortedNeighbors(vRemove, scratchNbrsA_);
        fillSortedNeighbors(vKeep, scratchNbrsB_);

        // Shared neighbors (excluding vKeep/vRemove from each other's sets)
        scratchSharedNbrs_.clear();
        for (auto v : scratchNbrsA_) {
            if (v != vKeep && std::binary_search(scratchNbrsB_.begin(), scratchNbrsB_.end(), v)) {
                scratchSharedNbrs_.push_back(v);
            }
        }
        // scratchSharedNbrs_ is already sorted since scratchNbrsA_ is sorted

        // Expected shared: opposite vertices of the shared faces (at most 2 entries)
        std::array<std::size_t, 2> expectedShared{};
        std::size_t numExpected = 0;
        for (auto fi : sharedFaces) {
            for (auto vi : faces_[fi]) {
                if (vi != vRemove && vi != vKeep)
                    expectedShared[numExpected++] = vi;
            }
        }
        if (numExpected > 1 && expectedShared[0] > expectedShared[1])
            std::swap(expectedShared[0], expectedShared[1]);

        if (scratchSharedNbrs_.size() != numExpected)
            return std::nullopt;
        for (std::size_t i = 0; i < numExpected; ++i) {
            if (scratchSharedNbrs_[i] != expectedShared[i])
                return std::nullopt;
        }

        // Minimum angle threshold (radians) — reject collapses that would
        // create triangles with any angle below this.  Paper uses ~10°.
        constexpr T minAngle = PI<T> / T(18);  // 10°

        // Validate all faces that will exist around vKeep after collapse:
        // removeFaces (with vRemove→vKeep substitution) must not flip or
        // degenerate, and ALL surviving faces incident to vKeep must
        // maintain a minimum angle above the threshold.
        //
        // Collect the full set of post-collapse faces incident to vKeep.
        auto& postFaces = scratchPostFaces_;
        postFaces.clear();
        // Existing vKeep faces (excluding shared faces which will be removed)
        auto isInShared = [&](std::size_t fi) {
            return std::find(sharedFaces.begin(), sharedFaces.end(), fi) != sharedFaces.end();
        };
        for (auto fi : vertFaces_[vKeep]) {
            if (!faceAlive_[fi] || isInShared(fi)) {
                continue;
            }
            postFaces.push_back(faces_[fi]);
        }
        // removeFaces with vRemove→vKeep substitution
        for (auto fi : removeFaces) {
            std::array<std::size_t, 3> newTri = faces_[fi];
            for (auto& vi : newTri) {
                if (vi == vRemove) {
                    vi = vKeep;
                }
            }
            // Reject if substitution creates a degenerate face
            if (newTri[0] == newTri[1] || newTri[1] == newTri[2] || newTri[0] == newTri[2]) {
                return std::nullopt;
            }
            postFaces.push_back(newTri);

            // Also check for normal flip on the modified faces
            auto& op0 = positions_[faces_[fi][0]];
            auto& op1 = positions_[faces_[fi][1]];
            auto& op2 = positions_[faces_[fi][2]];
            auto oldNormal = cross(op1 - op0, op2 - op0);
            auto newNormal = cross(positions_[newTri[1]] - positions_[newTri[0]],
                                   positions_[newTri[2]] - positions_[newTri[0]]);
            if (dot(oldNormal, newNormal) < T(0)) {
                return std::nullopt;
            }
        }

        // Check all post-collapse faces for minimum angle
        for (auto& tri : postFaces) {
            auto& p0 = positions_[tri[0]];
            auto& p1 = positions_[tri[1]];
            auto& p2 = positions_[tri[2]];
            auto e01 = p1 - p0;
            auto e02 = p2 - p0;
            auto e12 = p2 - p1;
            auto l01 = norm(e01);
            auto l02 = norm(e02);
            auto l12 = norm(e12);
            if (l01 == T(0) || l02 == T(0) || l12 == T(0)) {
                return std::nullopt;
            }
            // Clamp acos argument to [-1,1] for numerical safety
            auto clampedAngle = [](T cosVal) -> T {
                return std::acos(std::max(T(-1), std::min(T(1), cosVal)));
            };
            T a0 = clampedAngle(dot(e01, e02) / (l01 * l02));
            T a1 = clampedAngle(dot(p0 - p1, e12) / (l01 * l12));
            T a2 = PI<T> - a0 - a1;
            if (a0 < minAngle || a1 < minAngle || a2 < minAngle) {
                return std::nullopt;
            }
        }

        // Build collapse record with barycentric coordinates
        // After collapse, face (vRemove, vA, vB) becomes (vKeep, vA, vB).
        // Store bary coords of vRemoved's position in the post-collapse triangle.
        CollapseRecord<T> record;
        record.vRemoved = vRemove;
        record.vKept = vKeep;

        if (!removeFaces.empty()) {
            // Use the first surviving face; its post-collapse vertices are
            // (vKeep, vA, vB) where vA and vB are the non-vRemove vertices.
            auto fi = removeFaces[0];
            std::array<std::size_t, 3> postTri;
            postTri[0] = vKeep;
            std::size_t slot = 1;
            for (auto vi : faces_[fi]) {
                if (vi != vRemove) {
                    postTri[slot++] = vi;
                }
            }
            record.containingTri = postTri;
            record.bary = computeBarycentric_(positions_[postTri[0]], positions_[postTri[1]],
                                              positions_[postTri[2]], positions_[vRemove]);
        } else {
            // Edge case: all faces are shared — vertex collapses directly onto vKeep
            record.containingTri = {vKeep, vKeep, vKeep};
            record.bary = {T(1), T(0), T(0)};
        }

        // Execute collapse
        // Kill shared faces
        for (auto fi : sharedFaces) {
            faceAlive_[fi] = false;
            numAliveFaces_--;
        }

        // Update removeFaces: replace vRemove with vKeep
        for (auto fi : removeFaces) {
            for (auto& vi : faces_[fi]) {
                if (vi == vRemove) {
                    vi = vKeep;
                }
            }
            vertFaces_[vKeep].push_back(fi);
        }

        // Mark vRemove as dead
        alive_[vRemove] = false;
        numAliveVerts_--;

        // Propagate boundary status: if vRemove was on the boundary,
        // vKeep inherits it (it now sits on the mesh boundary).
        if (isBoundary_[vRemove]) {
            isBoundary_[vKeep] = true;
        }

        // Merge quadrics
        quadrics_[vKeep] += quadrics_[vRemove];

        // Compact dead face indices from vKeep's adjacency list
        auto& vkFaces = vertFaces_[vKeep];
        vkFaces.erase(std::remove_if(vkFaces.begin(), vkFaces.end(),
                                     [this](std::size_t fi) { return !faceAlive_[fi]; }),
                      vkFaces.end());

        // Fill caller-owned post-collapse neighbor vector (reusing its
        // allocation across calls). Equivalent to vertexNeighbors(vKeep) but
        // without an additional heap allocation.
        if (outKeepNbrs != nullptr) {
            outKeepNbrs->clear();
            for (auto fi : vkFaces) {
                for (auto vi : faces_[fi]) {
                    if (vi != vKeep && alive_[vi]) {
                        outKeepNbrs->push_back(vi);
                    }
                }
            }
            std::sort(outKeepNbrs->begin(), outKeepNbrs->end());
            outKeepNbrs->erase(std::unique(outKeepNbrs->begin(), outKeepNbrs->end()),
                               outKeepNbrs->end());
        }

        return record;
    }

    /** Compute collapse cost for edge (v0 → v1): Q_merged evaluated at v1 */
    [[nodiscard]] auto collapseCost(std::size_t v0, std::size_t v1) const -> T
    {
        auto Q = quadrics_[v0] + quadrics_[v1];
        auto& p = positions_[v1];
        return Q.evaluate(p[0], p[1], p[2]);
    }

    /** Check if a vertex is alive */
    [[nodiscard]] auto isAlive(std::size_t v) const -> bool { return alive_[v]; }

    /** Check if a vertex is collapsible (not pinned, alive) */
    [[nodiscard]] auto isCollapsible(std::size_t v) const -> bool
    {
        return alive_[v] && !isPinned_[v];
    }

    /** Get edges incident to vertex v (pairs of (v, neighbor)) */
    [[nodiscard]] auto vertexNeighbors(std::size_t v) const -> std::vector<std::size_t>
    {
        std::vector<std::size_t> nbrs;
        for (auto fi : vertFaces_[v]) {
            if (!faceAlive_[fi])
                continue;
            for (auto vi : faces_[fi]) {
                if (vi != v && alive_[vi])
                    nbrs.push_back(vi);
            }
        }
        std::sort(nbrs.begin(), nbrs.end());
        nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
        return nbrs;
    }

    /** Take a snapshot of surviving vertices and faces for a hierarchy level */
    [[nodiscard]] auto snapshot() const -> HierarchyLevel<T>
    {
        HierarchyLevel<T> level;
        level.originalToLocal.assign(alive_.size(), HierarchyLevel<T>::kAbsent);

        // Build mapping from original indices to level-local indices
        std::size_t localIdx = 0;
        for (std::size_t i = 0; i < alive_.size(); ++i) {
            if (alive_[i]) {
                level.originalToLocal[i] = localIdx;
                level.localToOriginal.push_back(i);
                level.positions.push_back(positions_[i]);
                localIdx++;
            }
        }

        // Remap faces
        for (std::size_t fi = 0; fi < faces_.size(); ++fi) {
            if (!faceAlive_[fi]) {
                continue;
            }
            std::array<std::size_t, 3> localTri;
            for (int j = 0; j < 3; ++j) {
                localTri[j] = level.originalToLocal[faces_[fi][j]];
            }
            level.faces.push_back(localTri);
        }

        return level;
    }

    /** Rebuild the edge list from alive faces and return it */
    auto rebuildAndGetEdges() -> const std::vector<std::pair<std::size_t, std::size_t>>&
    {
        buildEdges_();
        return edges_;
    }

private:
    /** Recompute per-vertex QEM quadrics from the live faces (Garland-Heckbert). */
    void computeQuadrics_()
    {
        for (auto& q : quadrics_) {
            q = Quadric<T>();
        }

        for (std::size_t fi = 0; fi < faces_.size(); ++fi) {
            if (!faceAlive_[fi]) {
                continue;
            }
            auto& tri = faces_[fi];
            auto& p0 = positions_[tri[0]];
            auto& p1 = positions_[tri[1]];
            auto& p2 = positions_[tri[2]];

            // Face plane: normal = (p1-p0) x (p2-p0), normalized
            auto e1 = p1 - p0;
            auto e2 = p2 - p0;
            auto n = cross(e1, e2);
            auto len = norm(n);
            if (len < std::numeric_limits<T>::epsilon()) {
                continue;
            }
            n /= len;

            T a = n[0], b = n[1], c = n[2];
            T d = -(a * p0[0] + b * p0[1] + c * p0[2]);

            Quadric<T> faceQ(a, b, c, d);
            quadrics_[tri[0]] += faceQ;
            quadrics_[tri[1]] += faceQ;
            quadrics_[tri[2]] += faceQ;
        }
    }

    /** Rebuild the unique-edge list from the live faces. */
    void buildEdges_()
    {
        edges_.clear();
        std::unordered_set<std::uint64_t> seen;
        auto edgeKey = [this](std::size_t a, std::size_t b) -> std::uint64_t {
            auto n = static_cast<std::uint64_t>(positions_.size());
            return static_cast<std::uint64_t>(std::min(a, b)) * n +
                   static_cast<std::uint64_t>(std::max(a, b));
        };

        for (std::size_t fi = 0; fi < faces_.size(); ++fi) {
            if (!faceAlive_[fi]) {
                continue;
            }
            auto& tri = faces_[fi];
            for (int j = 0; j < 3; ++j) {
                auto a = tri[j];
                auto b = tri[(j + 1) % 3];
                auto key = edgeKey(a, b);
                if (seen.insert(key).second) {
                    edges_.emplace_back(std::min(a, b), std::max(a, b));
                }
            }
        }
    }

    /** Compute barycentric coordinates of point p in triangle (a, b, c) */
    static auto computeBarycentric_(const Vec<T, 3>& a, const Vec<T, 3>& b, const Vec<T, 3>& c,
                                    const Vec<T, 3>& p) -> std::array<T, 3>
    {
        auto v0 = b - a;
        auto v1 = c - a;
        auto v2 = p - a;

        T d00 = dot(v0, v0);
        T d01 = dot(v0, v1);
        T d11 = dot(v1, v1);
        T d20 = dot(v2, v0);
        T d21 = dot(v2, v1);

        T denom = d00 * d11 - d01 * d01;
        if (std::abs(denom) < std::numeric_limits<T>::epsilon()) {
            return {T(1), T(0), T(0)};
        }

        T v = (d11 * d20 - d01 * d21) / denom;
        T w = (d00 * d21 - d01 * d20) / denom;

        // Clamp to avoid extreme extrapolation from far-away collapses
        v = std::max(T(-0.5), std::min(T(1.5), v));
        w = std::max(T(-0.5), std::min(T(1.5), w));
        T u = T(1) - v - w;

        return {u, v, w};
    }

    /** Per-vertex 3D positions, indexed by original vertex index. */
    std::vector<Vec<T, 3>> positions_;
    /** Per-vertex alive flag; collapses set entries to false. */
    std::vector<bool> alive_;
    /** Per-vertex boundary flag, copied from the source mesh. */
    std::vector<bool> isBoundary_;
    /** Per-vertex pinned flag; pinned vertices cannot be removed by collapse. */
    std::vector<bool> isPinned_;
    /** Per-vertex QEM quadrics accumulated from incident face planes. */
    std::vector<Quadric<T>> quadrics_;
    /** Face connectivity: each face is three vertex indices into `positions_`. */
    std::vector<std::array<std::size_t, 3>> faces_;
    /** Per-face alive flag; collapses retire faces by setting entries to false. */
    std::vector<bool> faceAlive_;
    /** Per-vertex incident face list. */
    std::vector<std::vector<std::size_t>> vertFaces_;
    /** Count of alive vertices (kept current across collapses). */
    std::size_t numAliveVerts_{0};
    /** Count of alive faces (kept current across collapses). */
    std::size_t numAliveFaces_{0};
    /** Unique-edge list, rebuilt on demand by `buildEdges_()`. */
    std::vector<std::pair<std::size_t, std::size_t>> edges_;

    /** Scratch buffer: neighbor indices shared between two endpoints. */
    std::vector<std::size_t> scratchShared_;
    /** Scratch buffer: faces to retire on a collapse. */
    std::vector<std::size_t> scratchRemove_;
    /** Scratch buffer: post-collapse face rewrites. */
    std::vector<std::array<std::size_t, 3>> scratchPostFaces_;
    /** Scratch buffer: neighbor list of one endpoint. */
    std::vector<std::size_t> scratchNbrsA_;
    /** Scratch buffer: neighbor list of the other endpoint. */
    std::vector<std::size_t> scratchNbrsB_;
    /** Scratch buffer: deduplicated union of `scratchNbrsA_` and `scratchNbrsB_`. */
    std::vector<std::size_t> scratchSharedNbrs_;
};

/**
 * @brief Build a mesh hierarchy by greedy QEM decimation
 *
 * Returns a vector of HierarchyLevel from finest to coarsest, plus
 * the collapse records needed for prolongation (ordered from finest to coarsest).
 */
template <typename T, class MeshPtr>
auto buildHierarchy(const MeshPtr& mesh, std::size_t pin0, std::size_t pin1, std::size_t levelRatio,
                    std::size_t minCoarseVerts)
    -> std::pair<std::vector<HierarchyLevel<T>>, std::vector<std::vector<CollapseRecord<T>>>>
{
    DecimationMesh<T> dmesh;
    dmesh.build(mesh, pin0, pin1);

    // Finest level snapshot
    std::vector<HierarchyLevel<T>> levels;
    levels.push_back(dmesh.snapshot());

    std::vector<std::vector<CollapseRecord<T>>> collapsesByLevel;

    auto targetVerts = dmesh.numAliveVerts();
    if (targetVerts <= minCoarseVerts) {
        // Mesh is already small enough — single level
        return {levels, collapsesByLevel};
    }

    // Reused across collapses to avoid per-collapse heap allocation; populated
    // by tryCollapse with vKeep's post-collapse neighbors.
    std::vector<std::size_t> nbrs;

    while (targetVerts > minCoarseVerts) {
        auto nextTarget = std::max(targetVerts / levelRatio, minCoarseVerts);
        std::vector<CollapseRecord<T>> levelCollapses;

        // Build priority queue of edge collapses
        using CostEdge = std::pair<T, std::pair<std::size_t, std::size_t>>;
        std::priority_queue<CostEdge, std::vector<CostEdge>, std::greater<CostEdge>> pq;

        const auto& edges = dmesh.rebuildAndGetEdges();
        for (auto& [a, b] : edges) {
            // Try collapsing the collapsible vertex towards the other
            if (dmesh.isCollapsible(a)) {
                pq.push({dmesh.collapseCost(a, b), {a, b}});
            }
            if (dmesh.isCollapsible(b)) {
                pq.push({dmesh.collapseCost(b, a), {b, a}});
            }
        }

        while (dmesh.numAliveVerts() > nextTarget && !pq.empty()) {
            auto [cost, edge] = pq.top();
            pq.pop();

            auto [vRemove, vKeep] = edge;
            // Skip stale entries whose endpoints have already been collapsed
            if (!dmesh.isAlive(vRemove) || !dmesh.isAlive(vKeep)) {
                continue;
            }
            auto record = dmesh.tryCollapse(vRemove, vKeep, &nbrs);
            if (!record) {
                continue;
            }

            levelCollapses.push_back(*record);

            // Add new edges involving vKeep to the priority queue using the
            // post-collapse neighbor list populated by tryCollapse.
            for (auto nb : nbrs) {
                if (dmesh.isCollapsible(nb)) {
                    pq.push({dmesh.collapseCost(nb, vKeep), {nb, vKeep}});
                }
                if (dmesh.isCollapsible(vKeep)) {
                    pq.push({dmesh.collapseCost(vKeep, nb), {vKeep, nb}});
                }
            }
        }

        if (levelCollapses.empty()) {
            break;  // No more valid collapses possible
        }

        collapsesByLevel.push_back(std::move(levelCollapses));
        levels.push_back(dmesh.snapshot());
        targetVerts = dmesh.numAliveVerts();
    }

    return {levels, collapsesByLevel};
}

/**
 * @brief Build a HalfEdgeMesh from a hierarchy level
 */
template <typename T>
auto buildLevelMesh(const HierarchyLevel<T>& level) -> typename HalfEdgeMesh<T>::Pointer
{
    auto mesh = HalfEdgeMesh<T>::New();

    for (const auto& pos : level.positions) {
        mesh->insert_vertex(pos[0], pos[1], pos[2]);
    }

    // level.faces is already `vector<array<size_t,3>>`; insert_faces is generic
    // over containers-of-iterables so we can pass it directly instead of
    // copying every face into a `vector<vector<size_t>>` (one heap allocation
    // per face).
    mesh->insert_faces(level.faces);

    return mesh;
}

/** Sentinel value marking a vertex slot in a UV vector that has not yet been
 *  solved/prolongated. NaN propagates if read by accident, making mistakes
 *  loud. */
template <typename T>
constexpr auto kUnsetUV =
    std::array<T, 2>{std::numeric_limits<T>::quiet_NaN(), std::numeric_limits<T>::quiet_NaN()};

/**
 * @brief Prolongate UV coordinates from a coarser level to a finer level
 *
 * Surviving vertices get their UVs directly; removed vertices get UVs
 * via barycentric interpolation in their containing post-collapse triangle.
 *
 * @param coarseUVs UV coordinates indexed by original vertex index
 *                  (unset slots contain `kUnsetUV<T>` / NaN).
 * @param collapses Collapse records for this level transition (finest-to-coarsest order)
 * @return UV vector indexed by original vertex index (includes all finer-level vertices
 *         touched by collapses; other slots remain `kUnsetUV<T>`).
 */
template <typename T>
auto prolongateUVs(std::vector<std::array<T, 2>> coarseUVs,
                   const std::vector<CollapseRecord<T>>& collapses) -> std::vector<std::array<T, 2>>
{
    // Start with the coarse-level UVs and fill in vRemoved slots in reverse.
    auto& fineUVs = coarseUVs;

    // Undo collapses in reverse order (coarsest collapse first was last applied)
    for (auto it = collapses.rbegin(); it != collapses.rend(); ++it) {
        auto& rec = *it;
        auto& tri = rec.containingTri;

        // All three containing-tri vertices should have UVs by now
        const auto& uv0 = fineUVs[tri[0]];
        const auto& uv1 = fineUVs[tri[1]];
        const auto& uv2 = fineUVs[tri[2]];

        std::array<T, 2> newUV;
        newUV[0] = rec.bary[0] * uv0[0] + rec.bary[1] * uv1[0] + rec.bary[2] * uv2[0];
        newUV[1] = rec.bary[0] * uv0[1] + rec.bary[1] * uv1[1] + rec.bary[2] * uv2[1];
        fineUVs[rec.vRemoved] = newUV;
    }

    return fineUVs;
}

/**
 * @brief Solve the LSCM system at one hierarchy level
 *
 * Builds the Lévy et al. Eq. 10 LSCM system on the given level mesh with the
 * given pin vertices. If an initial guess is provided, uses solveWithGuess.
 *
 * @param origVertCount Original (finest-level) vertex count; sizes the output
 *                      UV vector so it can be indexed by original vertex idx.
 * @return UV coordinates indexed by original vertex index. Slots for vertices
 *         not present at this level remain `kUnsetUV<T>` (NaN).
 */
template <typename T, class SolverType>
auto solveLSCMLevel(const typename HalfEdgeMesh<T>::Pointer& levelMesh,
                    const detail::hlscm::HierarchyLevel<T>& level, std::size_t origPin0,
                    std::size_t origPin1, std::size_t origVertCount,
                    const std::vector<std::array<T, 2>>* initialGuess)
    -> std::vector<std::array<T, 2>>
{
    using SparseMatrix = Eigen::SparseMatrix<T>;
    using DenseMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using Mesh = HalfEdgeMesh<T>;

    auto numVerts = levelMesh->num_vertices();
    constexpr std::size_t numFixed = 2;
    auto numFree = numVerts - numFixed;

    // Map original pin indices to level-local indices
    auto localPin0 = level.originalToLocal[origPin0];
    auto localPin1 = level.originalToLocal[origPin1];
    auto p0 = levelMesh->vertex(localPin0);
    auto p1 = levelMesh->vertex(localPin1);

    // Pin placement + LSCM system assembly (shared with AngleBasedLSCM)
    auto parts = detail::lscm::buildSystem<T, Mesh>(levelMesh, p0, p1);
    auto& A = parts.A;
    auto& b = parts.b;
    auto& freeIdxTable = parts.freeIdxTable;

    // Build initial guess vector from prolongated UVs.
    // `initialGuess` is indexed by original vertex idx; entries with NaN
    // are vertices not yet solved at any coarser level.
    bool warmed = initialGuess && !initialGuess->empty();
    auto buildInitialGuess = [&]() -> DenseMatrix {
        DenseMatrix x0 = DenseMatrix::Zero(2 * numFree, 1);
        for (const auto& v : levelMesh->vertices()) {
            if (v == p0 || v == p1) {
                continue;
            }
            auto freeIdx = freeIdxTable.at(v->idx);
            auto origIdx = level.localToOriginal[v->idx];
            const auto& guess = (*initialGuess)[origIdx];
            if (!std::isnan(guess[0])) {
                x0(2 * freeIdx, 0) = guess[0];
                x0(2 * freeIdx + 1, 0) = guess[1];
            }
        }
        return x0;
    };

    // Convergence tolerance for iterative solvers.  Eigen defaults to
    // machine epsilon (~2e-16 for double) which is far tighter than
    // needed for UV parameterization and prevents the warm-start from
    // reducing iteration count.  1e-8 gives ~8 digits of relative
    // residual precision — more than sufficient for texturing.
    constexpr T kTolerance = T(1e-8);

    // Solve
    DenseMatrix x;
    if constexpr (detail::is_instance_of_v<SolverType, Eigen::LeastSquaresConjugateGradient>) {
        // LSCG operates on the rectangular system A directly (avoids squaring the condition number)
        SolverType solver(A);
        solver.setTolerance(kTolerance);
        DenseMatrix bDense = DenseMatrix(b);
        if (warmed) {
            x = solver.solveWithGuess(bDense, buildInitialGuess());
        } else {
            x = solver.solve(bDense);
        }
        if (solver.info() == Eigen::ComputationInfo::NumericalIssue ||
            solver.info() == Eigen::ComputationInfo::InvalidInput ||
            solver.info() == Eigen::ComputationInfo::NoConvergence) {
            throw SolverException("HLSCM: LSCG solve failed at hierarchy level");
        }
    } else if constexpr (std::is_base_of_v<Eigen::IterativeSolverBase<SolverType>, SolverType>) {
        // CG and other iterative solvers on the square SPD system AtA.
        SparseMatrix AtA = A.transpose() * A;
        AtA.makeCompressed();
        DenseMatrix Atb = DenseMatrix(A.transpose() * b);
        SolverType solver(AtA);
        solver.setTolerance(kTolerance);
        if (warmed) {
            x = solver.solveWithGuess(Atb, buildInitialGuess());
        } else {
            x = solver.solve(Atb);
        }
        if (solver.info() == Eigen::ComputationInfo::NumericalIssue ||
            solver.info() == Eigen::ComputationInfo::InvalidInput ||
            solver.info() == Eigen::ComputationInfo::NoConvergence) {
            throw SolverException("HLSCM: iterative solve failed at hierarchy level");
        }
    } else {
        // Direct solver: decompose AtA and solve. No warm-start benefit.
        SparseMatrix AtA = A.transpose() * A;
        AtA.makeCompressed();
        DenseMatrix Atb = DenseMatrix(A.transpose() * b);
        SolverType solver;
        solver.compute(AtA);
        if (solver.info() != Eigen::ComputationInfo::Success) {
            throw SolverException("HLSCM: solver decomposition failed");
        }
        x = solver.solve(Atb);
    }

    // Build output UV vector indexed by original vertex idx. Vertices not
    // present at this level remain `kUnsetUV<T>`; they will be filled in by
    // prolongation when undoing collapses at finer levels.
    std::vector<std::array<T, 2>> uvs(origVertCount, kUnsetUV<T>);
    uvs[level.localToOriginal[p0->idx]] = {p0->pos[0], p0->pos[1]};
    uvs[level.localToOriginal[p1->idx]] = {p1->pos[0], p1->pos[1]};
    for (const auto& v : levelMesh->vertices()) {
        if (v == p0 || v == p1) {
            continue;
        }
        auto freeIdx = 2 * freeIdxTable.at(v->idx);
        auto origIdx = level.localToOriginal[v->idx];
        uvs[origIdx] = {x(freeIdx, 0), x(freeIdx + 1, 0)};
    }
    return uvs;
}

}  // namespace hlscm
}  // namespace detail

/**
 * @brief Compute parameterized mesh using Hierarchical LSCM
 *
 * Implements the HLSCM algorithm from Ray & Lévy, "Hierarchical Least Squares
 * Conformal Map" (2003) \cite ray2003hlscm. Uses cascadic multigrid to
 * accelerate LSCM
 * convergence: the mesh is decimated into a hierarchy, LSCM is solved on the
 * coarsest level, and the solution is prolongated and refined at each finer
 * level using conjugate gradient with the prolongated UVs as initial guess.
 *
 * For small meshes the hierarchy has a single level and HLSCM degrades
 * gracefully to a standard LSCM solve.
 *
 * @tparam T Floating-point type
 * @tparam MeshType HalfEdgeMesh type which implements the default mesh traits
 * @tparam Solver An Eigen iterative or direct solver. Iterative solvers
 *         (ConjugateGradient, LeastSquaresConjugateGradient) support warm-
 *         starting from the coarser-level solution; direct solvers ignore the
 *         initial guess. Defaults to
 *         `ConjugateGradient<SparseMatrix<T>, Lower|Upper>`, which solves the
 *         normal equations (AᵀA x = Aᵀb) of the overdetermined LSCM system.
 *         The `Lower|Upper` flag enables OpenMP-parallelized SpMV on the
 *         symmetric AᵀA matrix, giving the best multi-thread performance.
 *         LeastSquaresConjugateGradient is a valid alternative; it operates
 *         on the same normal equations internally but without the OpenMP
 *         benefit.
 *
 * @note **Solver default differs from AngleBasedLSCM.** AngleBasedLSCM
 *       defaults to SparseLU (a direct solver); HierarchicalLSCM defaults
 *       to ConjugateGradient so it can warm-start from the coarser-level
 *       solution. Using a direct solver via the `Solver` template parameter
 *       is valid but disables warm-starting — the initial guess is ignored.
 *
 * @note **Single-level fallback.** When the mesh is too small to decimate
 *       (all vertices are boundary, or minCoarseVertices is already reached),
 *       HLSCM falls back to a standard LSCM solve using *this class's*
 *       `Solver` template parameter, not AngleBasedLSCM's default (SparseLU).
 *       The result is numerically equivalent but may differ in convergence
 *       behavior from a plain AngleBasedLSCM call.
 */
template <typename T, class MeshType = HalfEdgeMesh<T>,
          class Solver =
              Eigen::ConjugateGradient<Eigen::SparseMatrix<T>, Eigen::Lower | Eigen::Upper>,
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

    /** @brief Set the vertex ratio between consecutive hierarchy levels (default: 10) */
    void setLevelRatio(std::size_t ratio)
    {
        if (ratio < 2) {
            throw std::invalid_argument("HierarchicalLSCM: levelRatio must be >= 2");
        }
        levelRatio_ = ratio;
    }

    /** @brief Set the minimum vertex count at the coarsest level (default: 100) */
    void setMinCoarseVertices(std::size_t count)
    {
        if (count < 3) {
            throw std::invalid_argument("HierarchicalLSCM: minCoarseVertices must be >= 3");
        }
        minCoarseVertices_ = count;
    }

    /**
     * @brief Compute parameterization using instance configuration
     *
     * Uses the pinned vertices, level ratio, and minimum coarse vertices
     * configured via `setPinnedVertices()`, `setLevelRatio()`, and
     * `setMinCoarseVertices()`. If no pins are set, selects them automatically
     * using the same boundary-walk logic as `Compute(mesh)`.
     *
     * @throws MeshException if pin selection fails (no boundary vertices)
     * @throws SolverException if any hierarchy level fails to solve
     */
    void compute(typename Mesh::Pointer& mesh) const
    {
        std::size_t p0, p1;
        if (pinnedVertices_) {
            p0 = pinnedVertices_->first;
            p1 = pinnedVertices_->second;
        } else {
            AutoSelectPins(mesh, p0, p1);
        }
        ComputeImpl(mesh, p0, p1, levelRatio_, minCoarseVertices_);
    }

    /**
     * @brief Compute with automatic pin selection
     *
     * Selects pins identically to AngleBasedLSCM::Compute().
     *
     * @throws MeshException if the mesh has no boundary vertices (pin selection
     *         fails) or the mesh is otherwise invalid
     * @throws SolverException if any hierarchy level fails to solve
     */
    static void Compute(typename Mesh::Pointer& mesh)
    {
        std::size_t p0, p1;
        AutoSelectPins(mesh, p0, p1);
        ComputeImpl(mesh, p0, p1);
    }

    /**
     * @brief Compute with explicit pinned vertex indices
     *
     * @throws SolverException if any hierarchy level fails to solve
     */
    static void Compute(typename Mesh::Pointer& mesh, std::size_t pin0Idx, std::size_t pin1Idx)
    {
        ComputeImpl(mesh, pin0Idx, pin1Idx);
    }

private:
    /** Select two pinned boundary vertices (same logic as AngleBasedLSCM) */
    static void AutoSelectPins(const typename Mesh::Pointer& mesh, std::size_t& p0, std::size_t& p1)
    {
        auto boundary = mesh->vertices_boundary();
        if (boundary.empty()) {
            throw MeshException("HierarchicalLSCM: mesh has no boundary vertices");
        }
        auto v0 = boundary.front();
        auto e = v0->edge;
        do {
            if (e->pair->is_boundary()) {
                break;
            }
            e = e->pair->next;
        } while (e != v0->edge);
        if (e == v0->edge && !e->pair->is_boundary()) {
            throw MeshException("Pinned vertex not on boundary");
        }
        p0 = v0->idx;
        p1 = e->next->vertex->idx;
    }

    /**
     * @brief Copy edge angles from the original mesh to a level mesh
     *
     * At the finest hierarchy level (k=0), the level mesh has the same
     * face/edge structure as the original mesh. If ABF was run beforehand,
     * we must use the ABF-optimized angles rather than recomputing from
     * geometry. Coarser levels always use geometry angles.
     */
    static void CopyAnglesFromOriginal(const typename Mesh::Pointer& original,
                                       const typename HalfEdgeMesh<T>::Pointer& levelMesh)
    {
        for (const auto& f : levelMesh->faces()) {
            auto origFace = original->face(f->idx);
            auto le = f->head;
            auto oe = origFace->head;
            for (int j = 0; j < 3; ++j) {
                le->alpha = oe->alpha;
                le = le->next;
                oe = oe->next;
            }
        }
    }

    /**
     * @brief Core hierarchical LSCM solve given resolved pin indices
     *
     * Builds the mesh hierarchy, solves LSCM at the coarsest level, then
     * prolongates and refines at each finer level.
     */
    static void ComputeImpl(typename Mesh::Pointer& mesh, std::size_t pin0Idx, std::size_t pin1Idx,
                            std::size_t levelRatio = 10, std::size_t minCoarseVerts = 100)
    {
        // Build mesh hierarchy
        auto [levels, collapsesByLevel] =
            detail::hlscm::buildHierarchy<T>(mesh, pin0Idx, pin1Idx, levelRatio, minCoarseVerts);

        if (levels.size() <= 1) {
            // Mesh too small for hierarchy — single-level LSCM solve
            // Use AngleBasedLSCM for exact equivalence on small meshes
            AngleBasedLSCM<T, MeshType, Solver>::Compute(mesh, pin0Idx, pin1Idx);
            return;
        }

        // UV vectors are sized to the original (finest-level) vertex count
        // so they can be indexed by `original vertex idx` at every level.
        const auto origVertCount = mesh->num_vertices();

        // Solve coarsest level (last in the array)
        auto coarsestIdx = levels.size() - 1;
        auto coarseMesh = detail::hlscm::buildLevelMesh<T>(levels[coarsestIdx]);
        ComputeMeshAngles(coarseMesh);
        auto uvs = detail::hlscm::solveLSCMLevel<T, Solver>(
            coarseMesh, levels[coarsestIdx], pin0Idx, pin1Idx, origVertCount, nullptr);

        // Prolongate and refine at each finer level
        for (std::size_t k = coarsestIdx; k-- > 0;) {
            // Prolongate UVs from level k+1 to level k
            uvs = detail::hlscm::prolongateUVs<T>(std::move(uvs), collapsesByLevel[k]);

            // Build level mesh
            auto levelMesh = detail::hlscm::buildLevelMesh<T>(levels[k]);

            if (k == 0) {
                // Finest level: use original mesh angles (may be ABF-optimized)
                CopyAnglesFromOriginal(mesh, levelMesh);
            } else {
                // Coarser levels: compute angles from 3D geometry
                ComputeMeshAngles(levelMesh);
            }

            // Solve with initial guess
            uvs = detail::hlscm::solveLSCMLevel<T, Solver>(levelMesh, levels[k], pin0Idx, pin1Idx,
                                                           origVertCount, &uvs);
        }

        // Transfer final UVs back to input mesh.
        // At the finest level every vertex has a UV (no NaN slots).
        for (const auto& v : mesh->vertices()) {
            const auto& uv = uvs[v->idx];
            v->pos[0] = uv[0];
            v->pos[1] = uv[1];
            v->pos[2] = T(0);
        }
    }

    /** Optional explicit pin pair */
    std::optional<std::pair<std::size_t, std::size_t>> pinnedVertices_;
    /** Ratio of vertices between consecutive hierarchy levels */
    std::size_t levelRatio_{10};
    /** Minimum vertex count for the coarsest hierarchy level */
    std::size_t minCoarseVertices_{100};
};

}  // namespace OpenABF


// #include "OpenABF/MeshIO.hpp"


#include <filesystem>
#include <fstream>

// #include "OpenABF/MeshIOFormats.hpp"


#include <charconv>
#include <filesystem>
#include <iostream>
#include <string_view>
#include <vector>

// #include "OpenABF/MeshIOUtils.hpp"


#include <algorithm>
#include <cctype>
#include <charconv>
#include <locale>
#include <string_view>
#include <vector>

namespace OpenABF::io_utils
{

/** @brief Compare two string_views, ignoring case */
static auto icase_compare(const std::string_view a, const std::string_view b) -> bool
{
    // not the same length
    if (a.length() != b.length()) {
        return false;
    }

    // iterate over the characters
    for (std::size_t i = 0; i < a.length(); ++i) {
        if (std::tolower(a[i]) != std::tolower(b[i])) {
            return false;
        }
    }

    // success
    return true;
}

/** @brief Left trim */
static auto trim_left(std::string_view s) -> std::string_view
{
    const auto& loc = std::locale();
    const auto start = std::find_if_not(std::begin(s), std::end(s),
                                        [&loc](auto ch) -> bool { return std::isspace(ch, loc); });
    s.remove_prefix(std::distance(std::begin(s), start));
    return s;
}

/** @brief Right trim */
static auto trim_right(std::string_view s) -> std::string_view
{
    const auto& loc = std::locale();
    const auto start = std::find_if_not(s.rbegin(), s.rend(), [&loc](auto ch) -> bool {
                           return std::isspace(ch, loc);
                       }).base();
    s.remove_suffix(std::distance(start, std::end(s)));
    return s;
}

/** @brief Trim from both ends */
static auto trim(std::string_view s) -> std::string_view
{
    s = trim_left(s);
    s = trim_right(s);
    return s;
}

/**
 * @brief Split a string by a delimiter
 *
 * When provided conflicting delimiters, the largest delimiter will take
 * precedence:
 *
 * ```{.cpp}
 * split("a->b->c", "-", "->");  // returns {"a", "b", "c"}
 * ```
 */
template <typename... Ds>
static auto split(std::string_view s, const Ds&... ds) -> std::vector<std::string_view>
{
    constexpr std::string_view DEFAULT_DELIM{" "};

    // Build delimiters list
    std::vector<std::string_view> delimiters;
    if (sizeof...(ds) > 0) {
        delimiters = {ds...};
    } else {
        delimiters.emplace_back(DEFAULT_DELIM);
    }

    // Get a list of all delimiter start pos and sizes
    std::vector<std::pair<std::string_view::size_type, std::string_view::size_type>> delimPos;
    for (const auto& delim : delimiters) {
        auto b = s.find(delim, 0);
        while (b != std::string_view::npos) {
            delimPos.emplace_back(b, delim.size());
            b = s.find(delim, b + delim.size());
        }
    }

    // Sort the delimiter start positions by first and largest
    std::sort(delimPos.begin(), delimPos.end(),
              [](const auto& l, const auto& r) { return l.second > r.second; });
    std::sort(delimPos.begin(), delimPos.end(),
              [](const auto& l, const auto& r) { return l.first < r.first; });

    // Split string
    std::vector<std::string_view> tokens;
    std::string_view::size_type begin{0};
    for (const auto& [end, size] : delimPos) {
        // ignore nested delimiters
        if (end < begin) {
            continue;
        }
        // get from begin to delim start
        if (auto t = s.substr(begin, end - begin); not t.empty()) {
            tokens.emplace_back(t);
        }
        begin = end + size;
    }
    if (auto t = s.substr(begin); not t.empty()) {
        tokens.emplace_back(t);
    }

    return tokens;
}

/**
 * @brief Convenience wrapper around std::to_chars for converting numerics to
 * std::string_view
 *
 * Useful during file writing operations when you're reusing a buffer, but
 * don't want to duplicate the error checking code of using std::to_chars.
 */
template <typename T>
auto to_string_view(const T& a, char* buf, const std::size_t& bufSize)
{
    auto res = std::to_chars(buf, buf + bufSize, a);
    if (res.ec != std::errc()) {
        throw std::runtime_error(std::make_error_code(res.ec).message());
    }
    return std::string_view(buf, res.ptr - buf);
}

/**
 * @brief Convert a string to a numeric type.
 *
 * A drop-in replacement for the `std:sto` family of functions which uses
 * `std::from_chars` for conversion. Like `std::sto`, throws exceptions when
 * conversion fails or if the converted value is out of range of the result
 * type.
 *
 * @throws std::invalid_argument If string cannot be converted to the result
 * type.
 * @throws std::result_out_of_range If converted value is out of range for the
 * result type.
 * @tparam T Requested numeric type
 * @tparam Args Parameter pack type
 * @param str Value to convert
 * @param args Extra parameters passed directly to `std::to_chars`
 * @return Converted value
 */
template <typename T, typename... Args>
auto to_numeric(const std::string_view str, Args... args) -> T
{
    T val;
    const auto* first = std::data(str);
    const auto* last = std::data(str) + std::size(str);
    auto [ptr, ec] = std::from_chars(first, last, val, args...);
    if (ec == std::errc::invalid_argument) {
        throw std::invalid_argument("Conversion could not be performed");
    }
    if (ec == std::errc::result_out_of_range) {
        throw std::out_of_range("Value out of range for the result type");
    }
    return val;
}

/**
 * @copybrief to_numeric
 *
 * Template specialization as fallback when the compiler does not support
 * `std::from_chars` for floating point types. Converts the input to a
 * `std::string` and passes to the appropriate `std::sto` function.
 */
template <>
inline auto to_numeric<float>(const std::string_view str) -> float
{
    return std::stof(std::string(str));
}

/** @copydoc to_numeric<float> */
template <>
inline auto to_numeric<double>(const std::string_view str) -> double
{
    return std::stod(std::string(str));
}

/** @copydoc to_numeric<float> */
template <>
inline auto to_numeric<long double>(const std::string_view str) -> long double
{
    return std::stold(std::string(str));
}
}  // namespace OpenABF::io_utils


namespace OpenABF::io_formats
{

/**
 * @brief Utility function for checking whether a given path matches one of the
 * accepted file extensions for a given format.
 */
template <typename PluginType>
static auto is_file_type(const std::filesystem::path& path)
{
    auto ext = path.extension().string();
    if (ext.empty()) {
        return false;
    }
    if (ext[0] == '.') {
        ext = ext.substr(1);
    }
    for (const auto& opt : PluginType::Extensions()) {
        if (io_utils::icase_compare(ext, opt)) {
            return true;
        }
    }
    return false;
}

/**
 * @brief Wavefront %OBJ file format
 *
 * Read support:
 *   - Vertices (v)
 *   - Faces (f) but only the first element (vertex ID) is stored
 *
 * Write support:
 *   - Vertices (v)
 *   - Vertex normals (vn)
 *   - Faces (f) but only v//vn constructions
 *
 * @see [Object Files (.obj) by Paul
 * Bourke](https://paulbourke.net/dataformats/obj/)
 */
struct OBJ {
    /** @brief List of recognized file format extensions */
    static auto Extensions() -> std::vector<std::string_view> { return {"obj"}; }

    /** Read the file stream into the provided object */
    template <typename MeshType>
    static auto Read(std::istream& is, MeshType& mesh)
    {
        using namespace io_utils;
        using T = typename MeshType::type;

        // Iterate the lines
        for (std::string line; std::getline(is, line);) {
            // Remove everything after a comment
            line = line.substr(0, line.find('#'));

            // Trim leading/trailing empty space
            auto line_view = trim(line);

            // Skip empty lines
            if (line_view.empty()) {
                continue;
            }

            // Split by part
            const auto parts = split(line_view);

            // Handle vertices
            if (parts[0] == "v") {
                std::vector<T> v;
                std::transform(parts.begin() + 1, parts.end(), std::back_inserter(v),
                               to_numeric<T>);
                mesh.insert_vertex(v);
            }

            // Handle faces (v attribute only)
            else if (parts[0] == "f") {
                std::vector<std::size_t> indices;
                std::transform(
                    parts.begin() + 1, parts.end(), std::back_inserter(indices),
                    [](const auto& p) { return to_numeric<std::size_t>(split(p, "/")[0]) - 1; });
                mesh.insert_face(indices);
            }
        }
        mesh.update_boundary();
    }

    /** Write the provided object to the given file stream */
    template <typename MeshType>
    static void Write(std::ostream& os, MeshType& mesh)
    {
        // Character buffer
        constexpr auto bufSize = 128;
        char buf[bufSize];

        // Write vertices
        for (std::size_t i = 0; i < mesh.num_vertices(); ++i) {
            const auto v = mesh.vertex(i);
            // write vertex position
            os << "v";
            for (const auto& a : v->pos) {
                auto res = std::to_chars(buf, buf + bufSize, a);
                if (res.ec != std::errc()) {
                    throw std::runtime_error(std::make_error_code(res.ec).message());
                }
                os << ' ' << std::string_view(buf, res.ptr - buf);
            }
            os << "\n";

            // write vertex normal
            os << "vn";
            for (const auto& a : v->normal()) {
                auto res = std::to_chars(buf, buf + bufSize, a);
                if (res.ec != std::errc()) {
                    throw std::runtime_error(std::make_error_code(res.ec).message());
                }
                os << ' ' << std::string_view(buf, res.ptr - buf);
            }
            os << "\n";
        }

        // Write faces
        for (std::size_t i = 0; i < mesh.num_faces(); ++i) {
            const auto f = mesh.face(i);
            os << "f";
            for (const auto& e : *f) {
                auto res = std::to_chars(buf, buf + bufSize, e->vertex->idx + 1);
                if (res.ec != std::errc()) {
                    throw std::runtime_error(std::make_error_code(res.ec).message());
                }
                // write vertex and normal IDs
                const auto id = std::string_view(buf, res.ptr - buf);
                os << ' ' << id << "//" << id;
            }
            os << "\n";
        }
    }
};

/**
 * @brief %PLY Polygon file format
 *
 * Read support:
 *   - Vertex properties:
 *     - (float) x, y, z
 *   - Face properties:
 *     - (property list uchar int) vertex_index
 *
 * Write support:
 *   - Formats: ASCII
 *   - Vertex properties:
 *     - (float) x, y, z, nx, ny, nz
 *   - Face properties:
 *     - (list uchar int) vertex_index
 *
 * @see [PLY - Polygon File Format by Paul
 * Bourke](https://paulbourke.net/dataformats/ply/)
 */
struct PLY {
    /** @brief List of recognized file format extensions */
    static auto Extensions() -> std::vector<std::string_view> { return {"ply"}; }

    /** Read the file stream into the provided object */
    template <typename MeshType>
    static auto Read(std::istream& is, MeshType& mesh)
    {
        using namespace io_utils;
        using T = typename MeshType::type;

        //// Parse header ////
        // Validate the type
        std::string line;
        std::getline(is, line);
        if (line != "ply") {
            throw std::runtime_error("File header does not begin with ply");
        }
        // Read the format line
        std::getline(is, line);
        const auto fmtParts = split(line);
        if (fmtParts[0] != "format") {
            throw std::runtime_error("File header missing format declaration");
        }
        if (fmtParts[1] != "ascii") {
            const auto fmt = std::string(fmtParts[1]) + " " + std::string(fmtParts[2]);
            throw std::runtime_error("Unsupported ply format: " + fmt);
        }

        // property = (label, type)
        struct Property {
            bool is_list{false};
            std::string list_count_type;
            std::string label;
            std::string type;
        };
        // element = (label, no. of elements, property list)
        struct Element {
            std::string label;
            std::uint32_t count{0};
            std::vector<Property> properties;
        };
        // list of elements
        std::vector<Element> elements;

        // Read the remaining header lines until end_header
        while (std::getline(is, line)) {
            // Trim leading/trailing empty space
            auto line_view = trim(line);

            // Skip empty lines
            if (line_view.empty()) {
                continue;
            }

            // Split by part
            const auto parts = split(line_view);

            // Handle comments (skip)
            if (parts[0] == "comment") {
                continue;
            }

            // Handle elements
            if (parts[0] == "element") {
                elements.push_back(
                    {.label = std::string(parts[1]), .count = to_numeric<std::uint32_t>(parts[2])});
            }

            // Handle properties for the most recent element
            else if (parts[0] == "property") {
                if (parts[1] == "list") {
                    elements.back().properties.push_back({.is_list = true,
                                                          .list_count_type = std::string(parts[2]),
                                                          .label = std::string(parts[4]),
                                                          .type = std::string(parts[3])});
                } else {
                    elements.back().properties.push_back(
                        {.label = std::string(parts[2]), .type = std::string(parts[1])});
                }
            }

            // Handle the end of the header
            else if (parts[0] == "end_header") {
                break;
            }
        }

        // Set up vertex map: v[n] -> property[m]
        // Probably unnecessary
        std::array<std::size_t, 3> vmap{};
        std::array<bool, 3> vmapFound{false, false, false};
        auto v_elem = std::find_if(elements.begin(), elements.end(),
                                   [](const auto& e) { return e.label == "vertex"; });
        if (v_elem == elements.end()) {
            throw std::runtime_error("Did not find vertex element");
        }
        for (auto i = 0; i < v_elem->properties.size(); ++i) {
            if (const auto& prop = v_elem->properties[i]; prop.label == "x") {
                vmap[0] = i;
                vmapFound[0] = true;
            } else if (prop.label == "y") {
                vmap[1] = i;
                vmapFound[1] = true;
            } else if (prop.label == "z") {
                vmap[2] = i;
                vmapFound[2] = true;
            }
        }
        if (!vmapFound[0] || !vmapFound[1] || !vmapFound[2]) {
            throw std::runtime_error("PLY vertex element missing required x/y/z properties");
        }

        // Iterate the lines of the body
        constexpr auto max_line = std::numeric_limits<std::streamsize>::max();
        for (const auto e : elements) {
            // Iterate the element lines
            for (auto i = 0; i < e.count; i++) {
                // parse vertex line
                if (e.label == "vertex") {
                    std::getline(is, line);
                    const auto line_view = trim(line);
                    const auto parts = split(line_view);
                    mesh.insert_vertex(to_numeric<T>(parts[vmap[0]]), to_numeric<T>(parts[vmap[1]]),
                                       to_numeric<T>(parts[vmap[2]]));
                }

                // parse face line
                else if (e.label == "face") {
                    std::getline(is, line);
                    const auto line_view = trim(line);
                    const auto parts = split(line_view);
                    if (parts[0] != "3") {
                        throw std::runtime_error("Unsupported number of vertices in face: " +
                                                 std::string(parts[0]));
                    }
                    mesh.insert_face(to_numeric<std::size_t>(parts[1]),
                                     to_numeric<std::size_t>(parts[2]),
                                     to_numeric<std::size_t>(parts[3]));
                }

                // ignore unrecognized element
                else {
                    is.ignore(max_line, is.widen('\n'));
                }
            }
        }
        mesh.update_boundary();
    }

    /** Write the provided object to the given file stream */
    template <typename MeshType>
    static void Write(std::ostream& os, MeshType& mesh)
    {
        using namespace io_utils;

        // Character buffer
        constexpr auto bufSize = 128;
        char buf[bufSize];

        // Write header
        os << "ply" << '\n';
        os << "format ascii 1.0" << '\n';
        os << "comment OpenABF PLY IO" << '\n';
        // Vertex element
        os << "element vertex ";
        os << to_string_view(mesh.num_vertices(), buf, bufSize) << '\n';
        os << "property float x" << '\n';
        os << "property float y" << '\n';
        os << "property float z" << '\n';
        os << "property float nx" << '\n';
        os << "property float ny" << '\n';
        os << "property float nz" << '\n';
        // Face element
        os << "element face ";
        os << to_string_view(mesh.num_faces(), buf, bufSize) << '\n';
        os << "property list uchar int vertex_indices" << '\n';
        os << "end_header" << '\n';

        // Write vertices
        for (std::size_t i = 0; i < mesh.num_vertices(); ++i) {
            const auto v = mesh.vertex(i);
            // write vertex position
            bool is_first{true};
            for (const auto& a : v->pos) {
                if (not is_first) {
                    os << ' ';
                }
                os << to_string_view(a, buf, bufSize);
                is_first = false;
            }
            for (const auto& a : v->normal()) {
                os << ' ' << to_string_view(a, buf, bufSize);
            }
            os << '\n';
        }

        // Write faces
        for (std::size_t i = 0; i < mesh.num_faces(); ++i) {
            // Only supports triangular faces
            os << '3';
            const auto f = mesh.face(i);
            for (const auto& e : *f) {
                os << ' ' << to_string_view(e->vertex->idx, buf, bufSize);
            }
            os << '\n';
        }
    }
};
}  // namespace OpenABF::io_formats

namespace OpenABF
{

/** @brief Load a HalfEdgeMesh from a file */
template <class MeshType>
auto ReadMesh(const std::filesystem::path& path)
{
    // Open the file
    std::ifstream file(path, std::ios::in);
    if (not file.is_open()) {
        throw std::runtime_error("Cannot open file for reading: " + path.string());
    }

    // Read the mesh
    auto result = MeshType::New();
    if (io_formats::is_file_type<io_formats::OBJ>(path)) {
        io_formats::OBJ::Read(file, *result);
    } else if (io_formats::is_file_type<io_formats::PLY>(path)) {
        io_formats::PLY::Read(file, *result);
    } else {
        throw std::runtime_error("Unsupported file type: " + path.extension().string());
    }

    return result;
}

/** @brief Write a HalfEdgeMesh to a file */
template <class MeshPtr>
void WriteMesh(const std::filesystem::path& path, const MeshPtr& mesh)
{
    // Open the file
    std::ofstream file(path, std::ios::out);
    if (not file.is_open()) {
        throw std::runtime_error("Cannot open file for writing: " + path.string());
    }

    // Write the mesh
    if (io_formats::is_file_type<io_formats::OBJ>(path)) {
        io_formats::OBJ::Write(file, *mesh);
    } else if (io_formats::is_file_type<io_formats::PLY>(path)) {
        io_formats::PLY::Write(file, *mesh);
    } else {
        throw std::runtime_error("Unsupported file type: " + path.extension().string());
    }

    // Close file
    file.flush();
    file.close();
    if (file.fail()) {
        throw std::runtime_error("Failed to write file: " + path.string());
    }
}

}  // namespace OpenABF

// clang-format on
