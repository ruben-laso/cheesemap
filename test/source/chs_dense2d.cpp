#include <cheesemap/cheesemap.hpp>

#include <gtest/gtest.h>

#include "mockup_data.hpp"

TEST(chs, chs_dense2d_build)
{
	auto points = chs::mockup::naive_points();

	EXPECT_NO_THROW(const auto map = chs::dense2d<chs::Point>(points, 1.0));
}

TEST(chs, chs_dense2d_build_random)
{
	static constexpr std::size_t num_builds = 100;
	static constexpr std::size_t num_points = 10'000;

	for ([[maybe_unused]] const auto _ : ranges::views::indices(num_builds))
	{
		const auto min_coord = get_random(-10.0, 10.0);
		const auto max_coord = get_random(90.0, 110.0);

		auto points = chs::mockup::random_points(num_points, { min_coord, min_coord, min_coord },
		                                         { max_coord, max_coord, max_coord });

		EXPECT_NO_THROW(const auto map = chs::dense2d<chs::Point>(points, 1.0));
	}
}

// TEST(chs, chs_dense2d_spherical_neighbours)
// {
// 	static constexpr std::size_t num_builds = 100;
// 	static constexpr std::size_t num_points = 10'000;
// 	static constexpr std::size_t num_neighbours = 10;
//
// 	for ([[maybe_unused]] const auto _ : ranges::views::indices(num_builds))
// 	{
// 		const auto min_coord = get_random(-10.0, 10.0);
// 		const auto max_coord = get_random(90.0, 110.0);
//
// 		auto points = chs::mockup::random_points(num_points, { min_coord, min_coord, min_coord },
// 		                                         { max_coord, max_coord, max_coord });
//
// 		const auto map = chs::dense2d<chs::Point>(points, 1.0);
//
// 		for (const auto & point : points)
// 		{
// 			const auto neighbours = map.spherical_neighbours(point, 1.0, num_neighbours);
// 			EXPECT_EQ(neighbours.size(), num_neighbours);
// 		}
// 	}
// }

auto main() -> int
{
	::testing::InitGoogleTest();
	return RUN_ALL_TESTS();
}