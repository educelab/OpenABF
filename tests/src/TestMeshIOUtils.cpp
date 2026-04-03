#include <gtest/gtest.h>
#include "OpenABF/OpenABF.hpp"

using namespace OpenABF::io_utils;

TEST(MeshIOUtils, ICaseCompare)
{
    EXPECT_TRUE(icase_compare("abc", "ABC"));
    EXPECT_TRUE(icase_compare("ABC", "abc"));
    EXPECT_TRUE(icase_compare("aBc", "AbC"));
    EXPECT_FALSE(icase_compare("abc", "abcd"));
    EXPECT_FALSE(icase_compare("abcd", "abc"));
    EXPECT_FALSE(icase_compare("abc", "abd"));
}

TEST(MeshIOUtils, TrimLeft)
{
    EXPECT_EQ(trim_left("  abc"), "abc");
    EXPECT_EQ(trim_left("abc  "), "abc  ");
    EXPECT_EQ(trim_left("  abc  "), "abc  ");
    EXPECT_EQ(trim_left(""), "");
    EXPECT_EQ(trim_left("   "), "");
}

TEST(MeshIOUtils, TrimRight)
{
    EXPECT_EQ(trim_right("  abc"), "  abc");
    EXPECT_EQ(trim_right("abc  "), "abc");
    EXPECT_EQ(trim_right("  abc  "), "  abc");
    EXPECT_EQ(trim_right(""), "");
    EXPECT_EQ(trim_right("   "), "");
}

TEST(MeshIOUtils, Trim)
{
    EXPECT_EQ(trim("  abc  "), "abc");
    EXPECT_EQ(trim("abc"), "abc");
    EXPECT_EQ(trim(""), "");
    EXPECT_EQ(trim("   "), "");
}

TEST(MeshIOUtils, Split)
{
    auto res = split("a b c");
    ASSERT_EQ(res.size(), 3);
    EXPECT_EQ(res[0], "a");
    EXPECT_EQ(res[1], "b");
    EXPECT_EQ(res[2], "c");

    res = split("a,b,c", ",");
    ASSERT_EQ(res.size(), 3);
    EXPECT_EQ(res[0], "a");
    EXPECT_EQ(res[1], "b");
    EXPECT_EQ(res[2], "c");

    res = split("a->b->c", "->");
    ASSERT_EQ(res.size(), 3);
    EXPECT_EQ(res[0], "a");
    EXPECT_EQ(res[1], "b");
    EXPECT_EQ(res[2], "c");

    res = split("a->b-c", "-", "->");
    ASSERT_EQ(res.size(), 3);
    EXPECT_EQ(res[0], "a");
    EXPECT_EQ(res[1], "b");
    EXPECT_EQ(res[2], "c");
}

TEST(MeshIOUtils, ToNumeric)
{
    EXPECT_EQ(to_numeric<int>("123"), 123);
    EXPECT_FLOAT_EQ(to_numeric<float>("123.45"), 123.45f);
    EXPECT_DOUBLE_EQ(to_numeric<double>("123.456"), 123.456);
    EXPECT_THROW(to_numeric<int>("abc"), std::invalid_argument);
    // std::out_of_range depends on T, but let's try a large number for int
    EXPECT_THROW(to_numeric<int>("99999999999999999999"), std::out_of_range);
}

TEST(MeshIOUtils, ToStringView)
{
    char buf[32];
    EXPECT_EQ(to_string_view(123, buf, 32), "123");
    EXPECT_EQ(to_string_view(123.45f, buf, 32), "123.45");
    EXPECT_THROW(to_string_view(123, buf, 1), std::runtime_error);
}
