/*
OpenABF
https://gitlab.com/educelab/OpenABF

Copyright 2021 EduceLab

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


#include <exception>
#include <string>

namespace OpenABF
{

/** @brief Solver exception */
class SolverException : std::exception
{
public:
    /** @brief Constructor with message */
    explicit SolverException(const char* msg) : msg_{msg} {}
    /** @brief Constructor with message */
    explicit SolverException(std::string msg) : msg_{std::move(msg)} {}
    /** @brief Get the exception message */
    const char* what() const noexcept override { return msg_.c_str(); }
private:
    /** Exception message */
    std::string msg_;
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
    return std::inner_product(
        std::begin(a), std::end(a), std::begin(b), Ret(0));
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
            return std::accumulate(
                std::begin(v), std::end(v), Ret(0),
                [](auto a, auto b) { return a + std::abs(b); });
        }
        case Norm::L2: {
            auto sum = std::accumulate(
                std::begin(v), std::end(v), Ret(0),
                [](auto a, auto b) { return a + (b * b); });
            return std::sqrt(sum);
        }
        case Norm::LInf: {
            return std::abs(*std::max_element(
                std::begin(v), std::end(v),
                [](auto a, auto b) { return std::abs(a) < std::abs(b); }));
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
template <
    typename T = float,
    typename T2,
    std::enable_if_t<std::is_floating_point<T>::value, bool> = true>
constexpr auto to_radians(T2 deg) -> T
{
    return deg * PI<T> / T(180);
}

/** @brief Convert radians to degrees */
template <
    typename T = float,
    typename T2,
    std::enable_if_t<std::is_floating_point<T>::value, bool> = true>
constexpr auto to_degrees(T2 rad) -> T
{
    return rad * T(180) / PI<T>;
}

}  // namespace OpenABF
// #include "OpenABF/Vector.hpp"


#include <array>

// #include "OpenABF/Math.hpp"


namespace OpenABF
{
namespace detail
{
/** @brief Helper type to perform parameter pack folding in C++11/14 */
struct ExpandType {
    /** Constructor */
    template <typename... T>
    explicit ExpandType(T&&...)
    {
    }
};
}  // namespace detail

/**
 * @brief N-dimensional vector class
 *
 * Essentially a wrapper around std::array that makes it more convenient for
 * vector math purposes.
 *
 * @tparam T Element type
 * @tparam Dims Number of elements
 */
template <
    typename T,
    std::size_t Dims,
    std::enable_if_t<std::is_arithmetic<T>::value, bool> = true>
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
#if __cplusplus >= 201703L
        // C++17 folding
        (val_[i++] = args, ...);
#elif __cplusplus > 201103L
        detail::ExpandType{0, ((val_[i++] = args), 0)...};
#endif
    }

    /** @brief Copy constructor */
    template <typename Vector>
    explicit Vec(const Vector& vec)
    {
        std::copy(val_.begin(), val_.end(), std::begin(vec));
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
    constexpr iterator rbegin() noexcept { return val_.rbegin(); }
    /** @brief Get an iterator to the first element of the vector */
    constexpr const_iterator rbegin() const noexcept { return val_.rbegin(); }
    /** @brief Get an iterator to the first element of the vector */
    constexpr const_iterator crbegin() const noexcept { return val_.crbegin(); }

    /**
     * @brief Get an iterator to one past the last element in the reverse vector
     */
    constexpr iterator rend() noexcept { return val_.rend(); }
    /**
     * @brief Get an iterator to one past the last element in the reverse vector
     */
    constexpr const_iterator rend() const noexcept { return val_.rend(); }
    /**
     * @brief Get an iterator to one past the last element in the reverse vector
     */
    constexpr const_iterator crend() const noexcept { return val_.crend(); }

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
    template <
        typename T2,
        std::enable_if_t<std::is_arithmetic<T2>::value, bool> = true>
    Vec& operator*=(const T2& b)
    {
        for (auto& v : val_) {
            v *= b;
        }
        return *this;
    }

    /** @brief Multiplication operator */
    template <class Vector>
    friend Vec operator*(Vec lhs, const Vector& rhs)
    {
        lhs *= rhs;
        return lhs;
    }

    /** @brief Division assignment operator */
    template <
        typename T2,
        std::enable_if_t<std::is_arithmetic<T2>::value, bool> = true>
    Vec& operator/=(const T2& b)
    {
        for (auto& v : val_) {
            v /= b;
        }
        return *this;
    }

    /** @brief Division operator */
    template <class Vector>
    friend Vec operator/(Vec lhs, const Vector& rhs)
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

}  // namespace OpenABF

/** Debug: Print a vector to a std::ostream */
template <typename T, std::size_t Dims>
std::ostream& operator<<(std::ostream& os, const OpenABF::Vec<T, Dims>& vec)
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

// #include "OpenABF/HalfEdgeMesh.hpp"


#include <algorithm>
#include <array>
#include <iterator>
#include <map>
#include <memory>
#include <vector>

// #include "OpenABF/Vector.hpp"


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

// #include "OpenABF/ABF.hpp"


#include <cmath>

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
struct ABFVertexTraits : public DefaultVertexTraits<T> {
    /** Lagrange Multiplier: Planarity constraint */
    T lambda_plan{0};
    /** Lagrange Multiplier: Reconstruction constraint */
    T lambda_len{1};
};

/** @brief ABF and ABFPlusPlus edge traits */
template <typename T>
struct ABFEdgeTraits : public DefaultEdgeTraits<T> {
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
struct ABFFaceTraits : public DefaultFaceTraits<T> {
    /** Lagrange Multiplier: Triangle validity constraint */
    T lambda_tri{0};
};
}  // namespace traits

namespace detail
{

/** @brief %ABF and ABF++ implementation details */
namespace ABF
{

/** @brief A HalfEdgeMesh with the %ABF traits */
template <typename T>
using Mesh = HalfEdgeMesh<
    T,
    3,
    traits::ABFVertexTraits<T>,
    traits::ABFEdgeTraits<T>,
    traits::ABFFaceTraits<T>>;

/** @brief Initialize the %ABF angles and weights from the edge alpha values */
template <typename T, class MeshPtr>
void InitializeAnglesAndWeights(MeshPtr& m)
{
    // Initialize and bound angle properties
    static constexpr auto MinAngle = PI<T> / T(180);
    static constexpr auto MaxAngle = PI<T> - MinAngle;
    for (auto& e : m->edges()) {
        e->alpha = e->beta = e->phi =
            std::min(std::max(e->alpha, MinAngle), MaxAngle);
        e->alpha_sin = std::sin(e->alpha);
        e->alpha_cos = std::cos(e->alpha);
        e->weight = T(1) / (e->phi * e->phi);
    }

    // Update weights for interior vertices
    for (auto& v : m->vertices_interior()) {
        auto wheel = v->wheel();
        auto angle_sum = std::accumulate(
            wheel.begin(), wheel.end(), T(0),
            [](auto a, auto b) { return a + b->beta; });
        for (auto& e : wheel) {
            e->phi *= 2 * PI<T> / angle_sum;
            e->weight = T(1) / (e->phi * e->phi);
        }
    }
}

/** @brief Compute ∇CTri w.r.t LambdaTri == CTri */
template <
    typename T,
    class FacePtr,
    std::enable_if_t<std::is_floating_point<T>::value, bool> = true>
auto TriGrad(const FacePtr& f) -> T
{
    T g = -PI<T>;
    for (const auto& e : *f) {
        g += e->alpha;
    }
    return g;
}

/** @brief Compute ∇CPlan w.r.t LambdaPlan == CPlan */
template <
    typename T,
    class VertPtr,
    std::enable_if_t<std::is_floating_point<T>::value, bool> = true>
auto PlanGrad(const VertPtr& v) -> T
{
    auto edges = v->wheel();
    T g = -2 * PI<T>;
    for (const auto& e : v->wheel()) {
        g += e->alpha;
    }
    return g;
}

/** @brief Compute ∇CLen w.r.t LambdaLen == CLen */
template <
    typename T,
    class VertPtr,
    std::enable_if_t<std::is_floating_point<T>::value, bool> = true>
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
template <
    typename T,
    class VertPtr,
    class EdgePtr,
    std::enable_if_t<std::is_floating_point<T>::value, bool> = true>
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
template <
    typename T,
    class EdgePtr,
    std::enable_if_t<std::is_floating_point<T>::value, bool> = true>
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
template <
    typename T,
    class MeshPtr,
    std::enable_if_t<std::is_floating_point<T>::value, bool> = true>
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
}  // namespace ABF
}  // namespace detail

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
template <
    typename T,
    class MeshType = detail::ABF::Mesh<T>,
    class Solver =
        Eigen::SparseLU<Eigen::SparseMatrix<T>, Eigen::COLAMDOrdering<int>>,
    std::enable_if_t<std::is_floating_point<T>::value, bool> = true>
class ABF
{
public:
    /** @brief Mesh type alias */
    using Mesh = MeshType;

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

    /** @brief Compute parameterized interior angles */
    void compute(typename Mesh::Pointer& mesh)
    {
        Compute(mesh, iters_, grad_, maxIters_);
    }

    /** @brief Compute parameterized interior angles */
    static void Compute(
        typename Mesh::Pointer& mesh,
        std::size_t& iters,
        T& gradient,
        std::size_t maxIters = 10)
    {
        using namespace detail::ABF;

        // Initialize angles and weights
        InitializeAnglesAndWeights<T>(mesh);

        // while ||∇F(x)|| > ε
        gradient = Gradient<T>(mesh);
        auto gradDelta = INF<T>;
        iters = 0;
        while (gradient > 0.001 and gradDelta > 0.001 and iters < maxIters) {
            // Typedefs
            using Triplet = Eigen::Triplet<T>;
            using SparseMatrix = Eigen::SparseMatrix<T>;
            using DenseVector = Eigen::Matrix<T, Eigen::Dynamic, 1>;

            // Helpful parameters
            auto vCnt = mesh->num_vertices();
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

            // vertex idx -> interior vertex idx permutation
            std::map<std::size_t, std::size_t> vIdx2vIntIdx;
            std::size_t newIdx{0};
            for (const auto& v : mesh->vertices_interior()) {
                vIdx2vIntIdx[v->idx] = newIdx++;
            }

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
                    triplets.emplace_back(row, e0->idx, 1);
                    triplets.emplace_back(e0->idx, row, 1);

                    // Jacobian of the CLen constraint
                    auto e1 = e0->next;
                    auto e2 = e1->next;
                    auto d1 = LenGrad<T>(v, e1);
                    auto d2 = LenGrad<T>(v, e2);
                    triplets.emplace_back(vIntCnt + row, e1->idx, d1);
                    triplets.emplace_back(vIntCnt + row, e2->idx, d2);
                    triplets.emplace_back(e1->idx, vIntCnt + row, d1);
                    triplets.emplace_back(e2->idx, vIntCnt + row, d2);
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
                throw SolverException(solver.lastErrorMessage());
            }
            DenseVector delta = solver.solve(b);
            if (solver.info() != Eigen::ComputationInfo::Success) {
                throw SolverException(solver.lastErrorMessage());
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
            for (auto& v : mesh->vertices_interior()) {
                auto intIdx = vIdx2vIntIdx.at(v->idx);
                v->lambda_plan += delta(idx + intIdx, 0);
                v->lambda_len += delta(idx + vIntCnt + intIdx, 0);
                idx++;
            }

            // Recalculate gradient for next iteration
            auto newGrad = detail::ABF::Gradient<T>(mesh);
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

protected:
    /** Gradient */
    T grad_{0};
    /** Number of executed iterations */
    std::size_t iters_{0};
    /** Max iterations */
    std::size_t maxIters_{10};
};

}  // namespace OpenABF
// #include "OpenABF/ABFPlusPlus.hpp"


#include <cmath>

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
template <
    typename T,
    class MeshType = detail::ABF::Mesh<T>,
    class Solver =
        Eigen::SparseLU<Eigen::SparseMatrix<T>, Eigen::COLAMDOrdering<int>>,
    std::enable_if_t<std::is_floating_point<T>::value, bool> = true>
class ABFPlusPlus
{
public:
    /** @brief Mesh type alias */
    using Mesh = MeshType;

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

    /** @brief Compute parameterized interior angles */
    void compute(typename Mesh::Pointer& mesh)
    {
        Compute(mesh, iters_, grad_, maxIters_);
    }

    /** @brief Compute parameterized interior angles */
    static void Compute(
        typename Mesh::Pointer& mesh,
        std::size_t& iters,
        T& gradient,
        std::size_t maxIters = 10)
    {
        using namespace detail::ABF;

        // Initialize angles and weights
        InitializeAnglesAndWeights<T>(mesh);

        // while ||∇F(x)|| > ε
        gradient = Gradient<T>(mesh);
        auto gradDelta = INF<T>;
        iters = 0;
        while (gradient > 0.001 and gradDelta > 0.001 and iters < maxIters) {
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

            // vertex idx -> interior vertex idx permutation
            std::map<std::size_t, std::size_t> vIdx2vIntIdx;
            std::size_t newIdx{0};
            for (const auto& v : mesh->vertices_interior()) {
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
            for (const auto& v : mesh->vertices_interior()) {
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
            for (auto& f : mesh->faces()) {
                f->lambda_tri += deltaLambda(f->idx, 0);
            }
            for (auto& v : mesh->vertices_interior()) {
                auto intIdx = vIdx2vIntIdx.at(v->idx);
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
};

}  // namespace OpenABF
// #include "OpenABF/AngleBasedLSCM.hpp"


#include <cmath>
#include <map>

#include <Eigen/SparseLU>

// #include "OpenABF/HalfEdgeMesh.hpp"

// #include "OpenABF/Math.hpp"


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
 * @tparam T Floating-point type
 * @tparam MeshType HalfEdgeMesh type which implements the default mesh traits
 * @tparam Solver A solver implementing the
 * [Eigen Sparse solver
 * concept](https://eigen.tuxfamily.org/dox-devel/group__TopicSparseSystems.html)
 * and templated on Eigen::SparseMatrix<T>
 */
template <
    typename T,
    class MeshType = HalfEdgeMesh<T>,
    class Solver =
        Eigen::SparseLU<Eigen::SparseMatrix<T>, Eigen::COLAMDOrdering<int>>,
    std::enable_if_t<std::is_floating_point<T>::value, bool> = true>
class AngleBasedLSCM
{
public:
    /** @brief Mesh type alias */
    using Mesh = MeshType;

    /** @brief Compute the parameterized mesh */
    void compute(typename Mesh::Pointer& mesh) const { Compute(mesh); }

    /** @brief Compute the parameterized mesh */
    static void Compute(typename Mesh::Pointer& mesh)
    {
        using Triplet = Eigen::Triplet<T>;
        using SparseMatrix = Eigen::SparseMatrix<T>;
        using DenseMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

        // Pinned vertex selection
        // Get the end points of a boundary edge
        auto p0 = mesh->vertices_boundary()[0];
        auto e = p0->edge;
        do {
            if (not e->pair) {
                break;
            }
            e = e->pair->next;
        } while (e != p0->edge);
        if (e == p0->edge and e->pair) {
            throw std::invalid_argument("Vertex not on boundary");
        }
        auto p1 = e->next->vertex;

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

            // If pin0 or pin1, put in fixedB matrix, else put in A
            auto row = 2 * f->idx;
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
// clang-format on