#include <cheesemap/maps/Sparse2D.hpp>

#include <gtest/gtest.h>

#include "mockup_data.hpp"

static constexpr auto CELL_SIZE = 5.0;

TEST(chs, chs_sparse2d_build)
{
	auto points = chs::mockup::naive_points();

	EXPECT_NO_THROW(const auto map = chs::Sparse2D<chs::Point>(points, CELL_SIZE));
}

TEST(chs, chs_sparse2d_build_random)
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

		EXPECT_NO_THROW(const auto map = chs::Sparse2D<chs::Point>(points, CELL_SIZE));
	}
}

TEST(chs, chs_sparse2d_spherical_search)
{
	static constexpr std::size_t num_builds = 20;
	static constexpr std::size_t num_points = 10'000;
	static constexpr std::size_t num_checks = 20;

	for ([[maybe_unused]] const auto _ : ranges::views::indices(num_builds))
	{
		const auto min_coord = get_random(chs::mockup::DEFAULT_MIN_COORD - chs::mockup::DEFAULT_MAP_VARIANCE,
		                                  chs::mockup::DEFAULT_MIN_COORD + chs::mockup::DEFAULT_MAP_VARIANCE);
		const auto max_coord = get_random(chs::mockup::DEFAULT_MAX_COORD - chs::mockup::DEFAULT_MAP_VARIANCE,
		                                  chs::mockup::DEFAULT_MAX_COORD + chs::mockup::DEFAULT_MAP_VARIANCE);

		auto points = chs::mockup::random_points(num_points, { min_coord, min_coord, min_coord },
		                                         { max_coord, max_coord, max_coord });

		auto map = chs::Sparse2D<chs::Point>(points, CELL_SIZE);

		for ([[maybe_unused]] const auto i : ranges::views::indices(num_checks))
		{
			const auto center = chs::mockup::random_point({ min_coord, min_coord, min_coord },
			                                              { max_coord, max_coord, max_coord });
			const auto radius = get_random(0.0, chs::mockup::DEFAULT_SEARCH_RADIUS);

			const auto sphere = chs::kernels::Sphere<3>(center, radius);

			const auto results_map = map.query(sphere);

			const auto results_truth =
			        ranges::views::filter(points, [&sphere](auto & p) { return sphere.is_inside(p); }) |
			        ranges::views::addressof | ranges::to<std::vector>();

			EXPECT_TRUE(are_the_same(results_map, results_truth));
		}
	}
}

TEST(chs, chs_sparse2d_spherical_search_filter)
{
	static constexpr std::size_t num_builds = 20;
	static constexpr std::size_t num_points = 10'000;
	static constexpr std::size_t num_checks = 20;

	for ([[maybe_unused]] const auto _ : ranges::views::indices(num_builds))
	{
		const auto min_coord = get_random(chs::mockup::DEFAULT_MIN_COORD - chs::mockup::DEFAULT_MAP_VARIANCE,
		                                  chs::mockup::DEFAULT_MIN_COORD + chs::mockup::DEFAULT_MAP_VARIANCE);
		const auto max_coord = get_random(chs::mockup::DEFAULT_MAX_COORD - chs::mockup::DEFAULT_MAP_VARIANCE,
		                                  chs::mockup::DEFAULT_MAX_COORD + chs::mockup::DEFAULT_MAP_VARIANCE);

		auto points = chs::mockup::random_points(num_points, { min_coord, min_coord, min_coord },
		                                         { max_coord, max_coord, max_coord });

		auto map = chs::Sparse2D<chs::Point>(points, CELL_SIZE);

		for ([[maybe_unused]] const auto i : ranges::views::indices(num_checks))
		{
			const auto center = chs::mockup::random_point({ min_coord, min_coord, min_coord },
			                                              { max_coord, max_coord, max_coord });
			const auto radius = get_random(0.0, chs::mockup::DEFAULT_SEARCH_RADIUS);

			const auto sphere = chs::kernels::Sphere<3>(center, radius);

			const auto filter = [](auto & p) { return p[2] > 50.0; };

			const auto results_map = map.query(sphere, filter);

			const auto results_truth =
			        ranges::views::filter(points, [&sphere](auto & p) { return sphere.is_inside(p); }) |
			        ranges::views::filter(filter) | ranges::views::addressof | ranges::to<std::vector>();

			EXPECT_TRUE(are_the_same(results_map, results_truth));
		}
	}
}


auto main() -> int
{
	::testing::InitGoogleTest();
	return RUN_ALL_TESTS();
}