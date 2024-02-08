#include <cheesemap/cheesemap.hpp>

#include <gtest/gtest.h>

TEST(chs, chs_exists)
{
	EXPECT_TRUE(true);
}

auto main() -> int
{
	::testing::InitGoogleTest();
	return RUN_ALL_TESTS();
}
