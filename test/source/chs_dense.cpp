#include <cheesemap/maps/Dense.hpp>

#include <gtest/gtest.h>

#include "generic_tests.hpp"
#include "mockup_data.hpp"

static constexpr auto CELL_SIZE = 5.0;

static constexpr auto FLAGS         = chs::flags::build::SHRINK_TO_FIT;
static constexpr auto REORDER_FLAGS = FLAGS | chs::flags::build::REORDER | chs::flags::build::PARALLEL;

const auto build_map         = [](auto & points) { return chs::Dense<chs::mockup::Point>(points, CELL_SIZE, FLAGS); };
const auto build_map_reorder = [](auto & points) {
	return chs::Dense<chs::mockup::Point>(points, CELL_SIZE, REORDER_FLAGS);
};

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
{ chs::test::build_random(build_map); }

TEST(chs, chs_dense_build_random_reorder)
{ chs::test::build_random(build_map_reorder); }

TEST(chs, chs_dense_sphere2d_search)
{ chs::test::kernel_search<chs::kernels::Sphere<chs::mockup::Point, 2>>(build_map); }

TEST(chs, chs_dense_sphere3d_search)
{ chs::test::kernel_search<chs::kernels::Sphere<chs::mockup::Point, 3>>(build_map); }

TEST(chs, chs_dense_sphere2d_search_filter)
{ chs::test::kernel_search_filter<chs::kernels::Sphere<chs::mockup::Point, 2>>(build_map); }

TEST(chs, chs_dense_sphere3d_search_filter)
{ chs::test::kernel_search_filter<chs::kernels::Sphere<chs::mockup::Point, 3>>(build_map); }

TEST(chs, chs_dense_sphere2d_search_reorder)
{ chs::test::kernel_search<chs::kernels::Sphere<chs::mockup::Point, 2>>(build_map_reorder); }

TEST(chs, chs_dense_sphere3d_search_reorder)
{ chs::test::kernel_search<chs::kernels::Sphere<chs::mockup::Point, 3>>(build_map_reorder); }

TEST(chs, chs_dense_sphere2d_search_filter_reorder)
{ chs::test::kernel_search_filter<chs::kernels::Sphere<chs::mockup::Point, 2>>(build_map_reorder); }

TEST(chs, chs_dense_sphere3d_search_filter_reorder)
{ chs::test::kernel_search_filter<chs::kernels::Sphere<chs::mockup::Point, 3>>(build_map_reorder); }

TEST(chs, chs_dense_cube_search)
{ chs::test::kernel_search<chs::kernels::Cube<chs::mockup::Point>>(build_map); }

TEST(chs, chs_dense_cube_search_filter)
{ chs::test::kernel_search_filter<chs::kernels::Cube<chs::mockup::Point>>(build_map); }

TEST(chs, chs_dense_cube_search_reorder)
{ chs::test::kernel_search<chs::kernels::Cube<chs::mockup::Point>>(build_map_reorder); }

TEST(chs, chs_dense_cube_search_filter_reorder)
{ chs::test::kernel_search_filter<chs::kernels::Cube<chs::mockup::Point>>(build_map_reorder); }

TEST(chs, chs_dense_knn)
{ chs::test::knn(build_map); }

TEST(chs, chs_dense_knn_reorder)
{ chs::test::knn(build_map_reorder); }

auto main() -> int
{
	::testing::InitGoogleTest();
	return RUN_ALL_TESTS();
}