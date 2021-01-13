#pragma once

#include <algorithm>
#include <array>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <vector>

namespace OpenABF
{

namespace detail
{
struct ExpandType {
    template <typename... T>
    explicit ExpandType(T&&...)
    {
    }
};
}  // namespace detail

template <typename T>
struct DefaultVertexProperties {
};

template <typename T>
struct DefaultEdgeProperties {
};

template <typename T>
struct DefaultFaceProperties {
};

template <
    typename T,
    size_t Dim = 3,
    typename VertexProperties = DefaultVertexProperties<T>,
    typename EdgeProperties = DefaultEdgeProperties<T>,
    typename FaceProperties = DefaultFaceProperties<T>>
class HalfEdgeMesh
{
private:
    struct Vertex;
    struct Edge;
    struct Face;

    using VertPtr = std::shared_ptr<Vertex>;
    using EdgePtr = std::shared_ptr<Edge>;
    using FacePtr = std::shared_ptr<Face>;

    struct Vertex : public VertexProperties {
        Vertex() = default;

        template <typename... Args>
        explicit Vertex(Args... args)
        {
            size_t i{0};
            detail::ExpandType{0, ((pos[i++] = args), 0)...};
        }

        template <typename... Args>
        static VertPtr New(Args&&... args)
        {
            return std::make_shared<Vertex>(std::forward<Args>(args)...);
        }

        auto edges() -> std::vector<EdgePtr>
        {
            std::vector<EdgePtr> ret;
            auto e = edge;
            do {
                if (not e->pair) {
                    throw std::runtime_error(
                        "Cannot enumerate edges for boundary_vertices vertex.");
                }
                ret.push_back(e);
                e = e->pair->next;
            } while (e != edge);
            return ret;
        }

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

        size_t idx{0};
        std::array<T, Dim> pos{};
        EdgePtr edge;
    };

    struct Edge : public EdgeProperties {
        template <typename... Args>
        static EdgePtr New(Args&&... args)
        {
            return std::make_shared<Edge>(std::forward<Args>(args)...);
        }

        auto is_boundary() const -> bool { return pair == nullptr; }

        EdgePtr pair;
        EdgePtr next;
        VertPtr vertex;
        FacePtr face;
    };

    struct Face : public FaceProperties {
        template <typename... Args>
        static FacePtr New(Args&&... args)
        {
            return std::make_shared<Face>(std::forward<Args>(args)...);
        }

        EdgePtr head;
    };

    std::vector<VertPtr> verts_;
    std::vector<FacePtr> faces_;
    std::multimap<size_t, EdgePtr> edges_;

public:
    HalfEdgeMesh() = default;

    template <typename... Args>
    auto insertVertex(Args... args) -> size_t
    {
        auto vert = Vertex::New(std::forward<Args>(args)...);
        vert->idx = verts_.size();
        verts_.push_back(vert);
        return vert->idx;
    }

    template <typename... Args>
    auto insertFace(Args... args) -> size_t
    {
        // Make a new face structure
        auto face = Face::New();

        // Iterate over the vertex indices
        size_t prevIdx;
        EdgePtr prevEdge;
        for (const auto& idx : {args...}) {
            // Make a new edge
            // TODO: Make sure edge doesn't already exist?
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

        faces_.emplace_back(face);
        return faces_.size() - 1;
    }

    auto vertices() -> std::vector<VertPtr> { return verts_; }

    auto edges() -> std::vector<EdgePtr> { return edges_; }

    auto faces() -> std::vector<FacePtr> { return faces_; }

    auto interior_vertices() -> std::vector<VertPtr>
    {
        std::vector<VertPtr> ret;
        std::copy_if(
            verts_.begin(), verts_.end(), std::back_inserter(ret),
            [](auto x) { return not x->is_boundary(); });
        return ret;
    }

    auto boundary_vertices() -> std::vector<VertPtr>
    {
        std::vector<VertPtr> ret;
        std::copy_if(
            verts_.begin(), verts_.end(), std::back_inserter(ret),
            [](auto x) { return x->is_boundary(); });
        return ret;
    }

private:
    // Find edge by start and end vertex indices
    auto find_edge_(size_t start, size_t end) -> EdgePtr
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