#include <sstream>
#include <gtest/gtest.h>
#include "OpenABF/OpenABF.hpp"
#include "Utils.hpp"

using namespace OpenABF;
using namespace OpenABF::tests;

using MeshType = HalfEdgeMesh<float>;

TEST(MeshIO, OBJ_ReadWrite)
{
    auto mesh = ConstructPyramid<MeshType>();

    std::stringstream ss;
    io_formats::OBJ::Write(ss, *mesh);

    auto mesh2 = MeshType::New();
    io_formats::OBJ::Read(ss, *mesh2);

    EXPECT_EQ(mesh->num_vertices(), mesh2->num_vertices());
    EXPECT_EQ(mesh->num_faces(), mesh2->num_faces());
}

TEST(MeshIO, PLY_ReadWrite)
{
    auto mesh = ConstructPyramid<MeshType>();

    std::stringstream ss;
    io_formats::PLY::Write(ss, *mesh);

    auto mesh2 = MeshType::New();
    io_formats::PLY::Read(ss, *mesh2);

    EXPECT_EQ(mesh->num_vertices(), mesh2->num_vertices());
    EXPECT_EQ(mesh->num_faces(), mesh2->num_faces());
}

TEST(MeshIO, OBJ_RoundTrip)
{
    auto mesh = ConstructPyramid<MeshType>();

    std::ostringstream oss;
    io_formats::OBJ::Write(oss, *mesh);

    auto mesh2 = MeshType::New();
    std::istringstream iss(oss.str());
    io_formats::OBJ::Read(iss, *mesh2);

    ASSERT_EQ(mesh->num_vertices(), mesh2->num_vertices());
    ASSERT_EQ(mesh->num_faces(), mesh2->num_faces());
    for (std::size_t v = 0; v < mesh->num_vertices(); ++v) {
        for (int i = 0; i < 3; ++i) {
            EXPECT_FLOAT_EQ(mesh->vertex(v)->pos[i], mesh2->vertex(v)->pos[i]);
        }
    }
}

TEST(MeshIO, PLY_RoundTrip)
{
    auto mesh = ConstructPyramid<MeshType>();

    std::ostringstream oss;
    io_formats::PLY::Write(oss, *mesh);

    auto mesh2 = MeshType::New();
    std::istringstream iss(oss.str());
    io_formats::PLY::Read(iss, *mesh2);

    ASSERT_EQ(mesh->num_vertices(), mesh2->num_vertices());
    ASSERT_EQ(mesh->num_faces(), mesh2->num_faces());
    for (std::size_t v = 0; v < mesh->num_vertices(); ++v) {
        for (int i = 0; i < 3; ++i) {
            EXPECT_FLOAT_EQ(mesh->vertex(v)->pos[i], mesh2->vertex(v)->pos[i]);
        }
    }
}

TEST(MeshIO, ReadWriteFile)
{
    auto mesh = ConstructPyramid<MeshType>();
    const std::string filename = "test_mesh.obj";

    EXPECT_NO_THROW(WriteMesh(filename, mesh));

    typename MeshType::Pointer mesh2;
    EXPECT_NO_THROW(mesh2 = ReadMesh<MeshType>(filename));

    EXPECT_EQ(mesh->num_vertices(), mesh2->num_vertices());
    EXPECT_EQ(mesh->num_faces(), mesh2->num_faces());

    std::filesystem::remove(filename);
}

TEST(MeshIO, UnsupportedFileType)
{
    auto mesh = ConstructPyramid<MeshType>();
    EXPECT_THROW(WriteMesh("test.unknown", mesh), std::runtime_error);
    EXPECT_THROW(ReadMesh<MeshType>("test.unknown"), std::runtime_error);
}

TEST(MeshIO, IsFileTypeNoExtension)
{
    EXPECT_FALSE(io_formats::is_file_type<io_formats::OBJ>("mesh"));
    EXPECT_FALSE(io_formats::is_file_type<io_formats::PLY>("mesh"));
}

TEST(MeshIO, PLY_MissingVertexProperty)
{
    // PLY header with x and y but missing z property
    std::string ply_data =
        "ply\n"
        "format ascii 1.0\n"
        "element vertex 1\n"
        "property float x\n"
        "property float y\n"
        "element face 0\n"
        "property list uchar int vertex_indices\n"
        "end_header\n"
        "1.0 2.0\n";

    std::istringstream iss(ply_data);
    auto mesh = MeshType::New();
    EXPECT_THROW(io_formats::PLY::Read(iss, *mesh), std::runtime_error);
}
