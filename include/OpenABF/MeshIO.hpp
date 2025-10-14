#pragma once

#include <filesystem>
#include <fstream>

#include "OpenABF/MeshIOFormats.hpp"

namespace OpenABF
{

/** @brief Load a HalfEdgeMesh from a file */
template <class MeshType>
auto ReadMesh(const std::filesystem::path& path)
{
    // Open the file
    std::ifstream file(path, std::ios::in);
    if (not file.is_open()) {
        throw std::runtime_error(
            "Cannot open file for reading: " + path.string());
    }

    // Read the mesh
    auto result = MeshType::New();
    if (io_formats::is_file_type<io_formats::OBJ>(path)) {
        io_formats::OBJ::Read(file, *result);
    } else if (io_formats::is_file_type<io_formats::PLY>(path)) {
        io_formats::PLY::Read(file, *result);
    } else {
        throw std::runtime_error(
            "Unsupported file type: " + path.extension().string());
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
        throw std::runtime_error(
            "Cannot open file for writing: " + path.string());
    }

    // Write the mesh
    if (io_formats::is_file_type<io_formats::OBJ>(path)) {
        io_formats::OBJ::Write(file, *mesh);
    } else if (io_formats::is_file_type<io_formats::PLY>(path)) {
        io_formats::PLY::Write(file, *mesh);
    } else {
        throw std::runtime_error(
            "Unsupported file type: " + path.extension().string());
    }

    // Close file
    file.flush();
    file.close();
    if (file.fail()) {
        throw std::runtime_error("Failed to write file: " + path.string());
    }
}

}  // namespace OpenABF
