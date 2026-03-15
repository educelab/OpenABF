#pragma once

#include <charconv>
#include <filesystem>
#include <iostream>
#include <string_view>
#include <vector>

#include "OpenABF/MeshIOUtils.hpp"

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