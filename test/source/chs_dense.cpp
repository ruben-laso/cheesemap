#include <cheesemap/maps/Dense.hpp>

#include <gtest/gtest.h>

#include "mockup_data.hpp"
#include "timing.hpp"

static constexpr auto CELL_SIZE = 5.0;

using map_type = chs::Dense<chs::Point, 2>;

TEST(chs, chs_dense_build)
{
	auto points = chs::mockup::naive_points();

	EXPECT_NO_THROW(const auto map = chs::Dense<chs::Point>(points, CELL_SIZE));
}

TEST(chs, chs_dense_build_random)
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

		EXPECT_NO_THROW(const auto map = map_type(points, CELL_SIZE));
	}
}

TEST(chs, chs_dense_spherical_search)
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

		auto map = map_type(points, CELL_SIZE);

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

TEST(chs, chs_dense_spherical_search_filter)
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

		auto map = map_type(points, CELL_SIZE);

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


TEST(chs, chs_dense_knn)
{
	static constexpr std::size_t num_builds = 20;
	static constexpr std::size_t num_points = 10'000;
	static constexpr std::size_t num_checks = 20;

	static constexpr std::size_t max_k = 2'000;

	double avg_secs_map   = 0.0;
	double avg_secs_truth = 0.0;

	for ([[maybe_unused]] const auto _ : ranges::views::indices(num_builds))
	{
		const auto min_coord = get_random(chs::mockup::DEFAULT_MIN_COORD - chs::mockup::DEFAULT_MAP_VARIANCE,
		                                  chs::mockup::DEFAULT_MIN_COORD + chs::mockup::DEFAULT_MAP_VARIANCE);
		const auto max_coord = get_random(chs::mockup::DEFAULT_MAX_COORD - chs::mockup::DEFAULT_MAP_VARIANCE,
		                                  chs::mockup::DEFAULT_MAX_COORD + chs::mockup::DEFAULT_MAP_VARIANCE);

		auto points = chs::mockup::random_points(num_points, { min_coord, min_coord, min_coord },
		                                         { max_coord, max_coord, max_coord });

		auto map = map_type(points, CELL_SIZE);

		for ([[maybe_unused]] const auto i : ranges::views::indices(num_checks))
		{
			const auto center = chs::mockup::random_point({ min_coord, min_coord, min_coord },
			                                              { max_coord, max_coord, max_coord });

			const auto distance = [](const auto & p, const auto & q) { return arma::norm(p - q); };

			const auto k = get_random({ 1 }, max_k);

			TIME_IT(const auto results_map = map.knn(k, center))
			avg_secs_map += chs::timing::last_seconds;

			auto get_truth = [&]() {
				using dist_point = std::pair<double, chs::Point *>;

				std::vector<dist_point> dist_points;
				dist_points.reserve(points.size());

				for (auto & p : points)
				{
					dist_points.emplace_back(distance(p, center), &p);
				}

				ranges::actions::sort(dist_points,
				                      [](const auto & a, const auto & b) { return a.first < b.first; });

				return dist_points | ranges::views::take(k) | ranges::to<std::vector>();
			};
			TIME_IT(const auto results_truth = get_truth())
			avg_secs_truth += chs::timing::last_seconds;

			// const auto same = are_the_same(results_map, results_truth);
			const auto same = are_the_same(results_map | ranges::views::values,
			                               results_truth | ranges::views::values);
			EXPECT_TRUE(same);
		}
	}

	std::cout << "Average seconds map: " << avg_secs_map / (num_builds * num_checks) << '\n';
	std::cout << "Average seconds truth: " << avg_secs_truth / (num_builds * num_checks) << '\n';
}

auto main() -> int
{
	::testing::InitGoogleTest();
	return RUN_ALL_TESTS();
}