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
