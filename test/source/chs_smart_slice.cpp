#include <cheesemap/maps/SmartSlice.hpp>

#include <gtest/gtest.h>

#include "mockup_data.hpp"

static constexpr auto CELL_SIZE = 5.0;

TEST(chs, chs_smart_slice_build)
{
	auto points = chs::mockup::naive_points();

	EXPECT_NO_THROW(const auto map = chs::slice::Smart<chs::Point>(points, CELL_SIZE));
}

TEST(chs, chs_smart_slice_build_random)
{
	static constexpr std::size_t num_builds = 100;
	static constexpr std::size_t num_points = 10'000;

	for ([[maybe_unused]] const auto _ : ranges::views::indices(num_builds))
	{
		const auto min_coord = get_random(chs::mockup::DEFAULT_MIN_COORD - chs::mockup::DEFAULT_MAP_VARIANCE,
		                                  chs::mockup::DEFAULT_MIN_COORD + chs::mockup::DEFAULT_MAP_VARIANCE);
		const auto max_coord = get_random(chs::mockup::DEFAULT_MAX_COORD - chs::mockup::DEFAULT_MAP_VARIANCE,
		                                  chs::mockup::DEFAULT_MAX_COORD + chs::mockup::DEFAULT_MAP_VARIANCE);

		auto points = chs::mockup::random_points(num_points, { min_coord, min_coord, min_coord },
		                                         { max_coord, max_coord, max_coord });

		EXPECT_NO_THROW(const auto map = chs::slice::Smart<chs::Point>(points, CELL_SIZE));
	}
}

TEST(chs, chs_smart_slice_sparsity)
{
	const auto cell_size = 5.0;
	chs::Point min{ 0.0, 0.0, 0.0 };
	chs::Point max{ 10.0, 10.0, 10.0 };
	chs::Box   box{ std::make_pair(min, max) };

	std::array sizes{ min[0] / cell_size, min[1] / cell_size };

	auto slice = chs::slice::Smart<chs::Point>(box, cell_size);

	EXPECT_TRUE(slice.sparse());
	EXPECT_DOUBLE_EQ(slice.density(), 0.0);

	// Insert point in cell[0, 0]
	auto point_00 = chs::mockup::random_point(min, { 5.0, 5.0, 5.0 });
	slice.add_point(point_00);
	EXPECT_TRUE(slice.sparse());
	EXPECT_DOUBLE_EQ(slice.density(), 0.25);

	// Insert point in cell[0, 1]
	auto point_01 = chs::mockup::random_point({ 0.0, 5.0, 0.0 }, { 5.0, 10.0, 5.0 });
	slice.add_point(point_01);
	EXPECT_TRUE(slice.sparse());
	EXPECT_DOUBLE_EQ(slice.density(), 0.5);

	// Insert point in cell[1, 0]
	auto point_10 = chs::mockup::random_point({ 5.0, 0.0, 0.0 }, { 10.0, 5.0, 5.0 });
	slice.add_point(point_10);
	EXPECT_TRUE(slice.sparse());
	EXPECT_DOUBLE_EQ(slice.density(), 0.75);

	// Insert point in cell[1, 1]
	auto point_11 = chs::mockup::random_point({ 5.0, 5.0, 0.0 }, { 10.0, 10.0, 5.0 });
	slice.add_point(point_11);
	EXPECT_FALSE(slice.sparse());
	EXPECT_DOUBLE_EQ(slice.density(), 1.0);
}

auto main() -> int
{
	::testing::InitGoogleTest();
	return RUN_ALL_TESTS();
}
