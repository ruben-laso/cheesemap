#include <cheesemap/maps/Dense.hpp>

#include <gtest/gtest.h>

#include "generic_tests.hpp"
#include "mockup_data.hpp"

static constexpr auto CELL_SIZE = 5.0;

auto build_map         = [](auto & points) { return chs::Dense<chs::Point>(points, CELL_SIZE, /* reorder = */ false); };
auto build_map_reorder = [](auto & points) { return chs::Dense<chs::Point>(points, CELL_SIZE, /* reorder = */ true); };

TEST(chs, chs_dense_build)
{
	auto points = chs::mockup::naive_points();

	EXPECT_NO_THROW(build_map(points));
}

TEST(chs, chs_dense_build_reorder)
{
	auto points = chs::mockup::naive_points();

	EXPECT_NO_THROW(build_map_reorder(points));
}

TEST(chs, chs_dense_build_random)
{
	chs::test::build_random(build_map);
}

TEST(chs, chs_dense_build_random_reorder)
{
	chs::test::build_random(build_map_reorder);
}

TEST(chs, chs_dense_spherical_search)
{
	chs::test::spherical_search(build_map);
}

TEST(chs, chs_dense_spherical_search_reorder)
{
	chs::test::spherical_search(build_map_reorder);
}

TEST(chs, chs_dense_spherical_search_filter)
{
	chs::test::spherical_search_filter(build_map);
}

TEST(chs, chs_dense_spherical_search_filter_reorder)
{
	chs::test::spherical_search_filter(build_map_reorder);
}

TEST(chs, chs_dense_knn)
{
	chs::test::knn(build_map);
}

TEST(chs, chs_dense_knn_reorder)
{
	chs::test::knn(build_map_reorder);
}

auto main() -> int
{
	::testing::InitGoogleTest();
	return RUN_ALL_TESTS();
}