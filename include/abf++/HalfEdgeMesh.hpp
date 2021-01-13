#pragma once

#include <algorithm>
#include <array>
#include <iostream>
#include <map>
#include <memory>
#include <vector>

namespace detail
{
struct ExpandType {
    template <typename... T>
    explicit ExpandType(T&&...)
    {
    }
};
}  // namespace detail

struct DefaultEdgeProperties {
};
struct DefaultVertexProperties {
};
struct DefaultFaceProperties {
};

template <
    typename T = float,
    size_t Dim = 3,
    typename VertexProperties = DefaultVertexProperties,
    typename EdgeProperties = DefaultEdgeProperties,
    typename FaceProperties = DefaultFaceProperties>
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

        bool is_boundary() const
        {
            if (edge) {
                return edge->pair == nullptr;
            }
            throw std::runtime_error(
                "Unreferenced vertex: " + std::to_string(idx));
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
    size_t insertVertex(Args... args)
    {
        auto vert = Vertex::New(std::forward<Args>(args)...);
        vert->idx = verts_.size();
        verts_.push_back(vert);
        return vert->idx;
    }

    template <typename... Args>
    size_t insertFace(Args... args)
    {
        // Make a new face structure
        auto face = Face::New();

        // Iterate over the vertex indices
        size_t prevIdx;
        EdgePtr prevEdge;
        for (const auto& idx : {args...}) {
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

private:
    // Find edge by start and end vertex indices
    EdgePtr find_edge_(size_t start, size_t end)
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