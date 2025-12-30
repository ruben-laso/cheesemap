#include "cheesemap/utils/morton.hpp"

#include <gtest/gtest.h>

TEST(Morton, Interleaving2D)
{
	// X = 5 (101), Y = 3 (011)
	// Expected: 1(X) 0(Y) 0(X) 1(Y) 1(X) 1(Y) -> 100111 (39 decimal)
	EXPECT_EQ(chs::morton::details::Dim2::morton(5, 3), 39);

    // X = 2 (10), Y = 1 (01)
    // Expected: 0(X) 0(Y) 1(X) 0(Y) 0(X) 1(Y) -> 001001 (9 decimal)
    EXPECT_EQ(chs::morton::details::Dim2::morton(2, 1), 9);
}

TEST(Morton, Boundary2D)
{
    // Check boundary (max 32-bit input)
    uint32_t max_val = 0xFFFFFFFF;
    uint64_t result  = chs::morton::details::Dim2::morton(max_val, max_val);
    EXPECT_EQ(result, 0xFFFFFFFFFFFFFFFF);
}

TEST(Morton, Interleaving3D)
{
	// X=1 (001), Y=1 (001), Z=1 (001)
	// Interleaved: X: _ _ 0  Y: _ 0 _  Z: 0 _ _ -> 000 000 111 (7 decimal)
	// Actually: (X<<2 | Y<<1 | Z) for each bit
    EXPECT_EQ(chs::morton::details::Dim3::morton(1, 1, 1), 7);

	// X=2(010), Y=0, Z=0 -> bits 5, 2, 1? No, just bit 3*1+2 = 5
	// X: 010 -> expand -> 100000 (32 decimal)
	EXPECT_EQ(chs::morton::details::Dim3::morton(2, 0, 0), 32);
}

TEST(Morton, Boundary3D)
{
    // Check boundary (max 21-bit input)
    uint32_t max_val = 0x1FFFFF; // 21 bits
    uint64_t result  = chs::morton::details::Dim3::morton(max_val, max_val, max_val);
    EXPECT_EQ(result, 0x7FFFFFFFFFFFFFFF);
}

TEST(MortonTest, Monotonicity2D) {
    auto m1 = chs::morton::details::Dim2::morton(10, 10);
    auto m2 = chs::morton::details::Dim2::morton(11, 10);
    auto m3 = chs::morton::details::Dim2::morton(11, 11);

    EXPECT_LT(m1, m2);
    EXPECT_LT(m2, m3);
}

TEST(MortonTest, Monotonicity3D) {
    auto m1 = chs::morton::details::Dim3::morton(10, 10, 10);
    auto m2 = chs::morton::details::Dim3::morton(11, 10, 10);
    auto m3 = chs::morton::details::Dim3::morton(11, 11, 10);
    auto m4 = chs::morton::details::Dim3::morton(11, 11, 11);

    EXPECT_LT(m1, m2);
    EXPECT_LT(m2, m3);
    EXPECT_LT(m3, m4);
}

TEST(Morton, Dim2)
{
	chs::Box   bbox({ 0.0, 0.0 }, { 10.0, 10.0 });
	chs::Point p1({ 0.1, 0.1 });
	chs::Point p2({ 1.0, 1.0 });
	chs::Point p3({ 0.5, 0.5 });
	chs::Point p4({ 5.0, 1.0 });
	chs::Point p5({ 10.0, 10.0 });

	auto code1 = chs::morton::details::Dim2::morton(p1, bbox);
	auto code2 = chs::morton::details::Dim2::morton(p2, bbox);
	auto code3 = chs::morton::details::Dim2::morton(p3, bbox);
	auto code4 = chs::morton::details::Dim2::morton(p4, bbox);
	auto code5 = chs::morton::details::Dim2::morton(p5, bbox);

	EXPECT_LT(code1, code3);
	EXPECT_LT(code3, code2);
	EXPECT_GT(code4, code2);
	EXPECT_GT(code5, code2);
}

TEST(Morton, Dim3)
{
	chs::Box   bbox({ 0.0, 0.0, 0.0 }, { 10.0, 10.0, 1.0 });
	chs::Point p1({ 0.1, 0.1, 0.1 });
	chs::Point p2({ 1.0, 1.0, 1.0 });
	chs::Point p3({ 0.5, 0.5, 0.5 });
	chs::Point p4({ 5.0, 1.0, 1.0 });
	chs::Point p5({ 10.0, 10.0, 1.0 });

	auto code1 = chs::morton::details::Dim3::morton(p1, bbox);
	auto code2 = chs::morton::details::Dim3::morton(p2, bbox);
	auto code3 = chs::morton::details::Dim3::morton(p3, bbox);
	auto code4 = chs::morton::details::Dim3::morton(p4, bbox);
	auto code5 = chs::morton::details::Dim3::morton(p5, bbox);

	EXPECT_LT(code1, code3);
	EXPECT_LT(code3, code2);
	EXPECT_GT(code4, code2);
	EXPECT_GT(code5, code2);
}

auto main() -> int
{
	::testing::InitGoogleTest();
	return RUN_ALL_TESTS();
}
