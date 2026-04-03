#include <cmath>
#include <sstream>
#include <gtest/gtest.h>

#include "OpenABF/OpenABF.hpp"

using namespace OpenABF;

TEST(Vec, Constructor)
{
    Vec3f a;
    EXPECT_EQ(a[0], 0.f);
    EXPECT_EQ(a[1], 0.f);
    EXPECT_EQ(a[2], 0.f);

    Vec3f b(1, 2, 3);
    Vec3f c(b);
    EXPECT_EQ(c, b);
}

TEST(Vec, Accessors)
{
    Vec3f a(1, 2, 3);
    EXPECT_EQ(a.at(0), 1.f);
    EXPECT_EQ(a.at(1), 2.f);
    EXPECT_EQ(a.at(2), 3.f);
    EXPECT_THROW(a.at(3), std::out_of_range);

    EXPECT_EQ(a.front(), 1.f);
    EXPECT_EQ(a.back(), 3.f);
    EXPECT_NE(a.data(), nullptr);
}

TEST(Vec, Modifiers)
{
    Vec3f a;
    a.fill(5.f);
    EXPECT_EQ(a, Vec3f(5, 5, 5));

    Vec3f b(1, 2, 3);
    a.swap(b);
    EXPECT_EQ(a, Vec3f(1, 2, 3));
    EXPECT_EQ(b, Vec3f(5, 5, 5));
}

TEST(Vec, Comparison)
{
    Vec3f a(1, 2, 3);
    Vec3f b(1, 2, 4);
    EXPECT_TRUE(a != b);
    EXPECT_FALSE(a == b);
}

TEST(Vec, Assignment)
{
    Vec3f a(1, 2, 3);
    Vec3f b;
    b = a;
    EXPECT_EQ(b, a);

    Vec3f c;
    c = {4.f, 5.f, 6.f};
    EXPECT_EQ(c, Vec3f(4, 5, 6));
}

TEST(Vec, Ostream)
{
    Vec3f a(1, 2, 3);
    std::stringstream ss;
    ss << a;
    EXPECT_EQ(ss.str(), "[1, 2, 3]");
}

TEST(Vec, Vec3d)
{
    Vec3d a(1.0, 2.0, 3.0);
    EXPECT_EQ(a.magnitude(), std::sqrt(14.0));
}

TEST(Vec, OperatorPlus)
{
    Vec3f a{1, 1, 1};
    Vec3f b{1, 1, 1};
    EXPECT_EQ(a + b, Vec3f(2, 2, 2));
    EXPECT_EQ(a, Vec3f(1, 1, 1));
    EXPECT_EQ(b, Vec3f(1, 1, 1));
    EXPECT_EQ(a += b, Vec3f(2, 2, 2));
    EXPECT_EQ(a, Vec3f(2, 2, 2));
}

TEST(Vec, OperatorMinus)
{
    Vec3f a{1, 1, 1};
    Vec3f b{1, 1, 1};
    EXPECT_EQ(a - b, Vec3f(0, 0, 0));
    EXPECT_EQ(a, Vec3f(1, 1, 1));
    EXPECT_EQ(b, Vec3f(1, 1, 1));
    EXPECT_EQ(a -= b, Vec3f(0, 0, 0));
    EXPECT_EQ(a, Vec3f(0, 0, 0));
}

TEST(Vec, OperatorMultiply)
{
    Vec3f a{1, 1, 1};
    EXPECT_EQ(a * 2, Vec3f(2, 2, 2));
    EXPECT_EQ(a, Vec3f(1, 1, 1));
    EXPECT_EQ(a *= 2, Vec3f(2, 2, 2));
    EXPECT_EQ(a, Vec3f(2, 2, 2));
}

TEST(Vec, OperatorDivide)
{
    Vec3f a{2, 2, 2};
    EXPECT_EQ(a / 2, Vec3f(1, 1, 1));
    EXPECT_EQ(a, Vec3f(2, 2, 2));
    EXPECT_EQ(a /= 2, Vec3f(1, 1, 1));
    EXPECT_EQ(a, Vec3f(1, 1, 1));
}

TEST(Vec, ScalarMultiplyDivide)
{
    EXPECT_EQ(Vec3f(2.f, 4.f, 6.f) * 2.f, Vec3f(4.f, 8.f, 12.f));
    EXPECT_EQ(Vec3f(4.f, 8.f, 12.f) / 2.f, Vec3f(2.f, 4.f, 6.f));
}

TEST(Vec, DotProduct)
{
    EXPECT_EQ(Vec3f(1, 0, 0).dot(Vec3f(0, 1, 0)), 0);
    EXPECT_EQ(Vec3f(1, 0, 0).dot(Vec3f(0, 0, 1)), 0);
    EXPECT_EQ(Vec3f(0, 1, 0).dot(Vec3f(0, 0, 1)), 0);

    EXPECT_EQ(Vec3f(1, 0, 0).dot(Vec3f(1, 0, 0)), 1);
    EXPECT_EQ(Vec3f(0, 1, 0).dot(Vec3f(0, 1, 0)), 1);
    EXPECT_EQ(Vec3f(0, 0, 1).dot(Vec3f(0, 0, 1)), 1);
}

TEST(Vec, CrossProduct)
{
    EXPECT_EQ(Vec3f(1, 0, 0).cross(Vec3f(1, 0, 0)), Vec3f(0, 0, 0));
    EXPECT_EQ(Vec3f(1, 0, 0).cross(Vec3f(0, 1, 0)), Vec3f(0, 0, 1));
    EXPECT_EQ(Vec3f(1, 0, 0).cross(Vec3f(0, 0, 1)), Vec3f(0, -1, 0));
}

TEST(Vec, Magnitude)
{
    EXPECT_EQ(Vec3f(1, 0, 0).magnitude(), 1.f);
    EXPECT_EQ(Vec3f(0, 2, 0).magnitude(), 2.f);
    EXPECT_EQ(Vec3f(0, 0, 3).magnitude(), 3.f);
}

TEST(Vec, UnitVector)
{
    Vec3f a{2, 0, 0};
    EXPECT_EQ(a.unit(), Vec3f(1, 0, 0));
    EXPECT_EQ(a, Vec3f(2, 0, 0));
}

TEST(Vec, ReverseIteration)
{
    Vec3f v{1.f, 2.f, 3.f};

    // rbegin/rend
    std::vector<float> result;
    for (auto it = v.rbegin(); it != v.rend(); ++it) {
        result.push_back(*it);
    }
    EXPECT_EQ(result, (std::vector<float>{3.f, 2.f, 1.f}));

    // crbegin/crend on a const Vec3f
    const Vec3f cv{1.f, 2.f, 3.f};
    std::vector<float> cresult;
    for (auto it = cv.crbegin(); it != cv.crend(); ++it) {
        cresult.push_back(*it);
    }
    EXPECT_EQ(cresult, (std::vector<float>{3.f, 2.f, 1.f}));
}
