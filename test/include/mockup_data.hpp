#pragma once

#include <armadillo>

#include <utility>
#include <vector>

#include "utils.hpp"

namespace chs::mockup
{
	using Point = arma::vec3;

	static constexpr double DEFAULT_MIN_COORD     = 0.0;
	static constexpr double DEFAULT_MAX_COORD     = 100.0;
	static constexpr double DEFAULT_MAP_VARIANCE  = 10.0;
	static constexpr double DEFAULT_SEARCH_RADIUS = 10.0;

	inline auto naive_points()
	{
		std::vector<Point> points = {
			{ 0.0, 0.0, 0.0 }, //
			{ 1.0, 1.0, 1.0 }, //
			{ 2.0, 2.0, 2.0 }, //
			{ 3.0, 3.0, 3.0 }, //
			{ 4.0, 4.0, 4.0 }, //
			{ 5.0, 5.0, 5.0 }, //
			{ 6.0, 6.0, 6.0 }, //
			{ 7.0, 7.0, 7.0 }, //
			{ 8.0, 8.0, 8.0 }, //
			{ 9.0, 9.0, 9.0 }  //
		};

		return points;
	}

	inline auto random_point(const Point & min, const Point & max)
	{
		Point point;

		for (const auto i : ranges::views::indices(3U))
		{
			point[i] = get_random(min[i], max[i]);
		}

		return point;
	}

	inline auto random_point()
	{
		const auto min = Point{ DEFAULT_MIN_COORD, DEFAULT_MIN_COORD, DEFAULT_MIN_COORD };
		const auto max = Point{ DEFAULT_MAX_COORD, DEFAULT_MAX_COORD, DEFAULT_MAX_COORD };

		return random_point(min, max);
	}

	inline auto random_points(const std::size_t n, const Point & min, const Point & max)
	{
		std::vector<Point> points(n);

		for (auto & point : points)
		{
			point = random_point(min, max);
		}

		return points;
	}

	inline auto random_points(const std::size_t n)
	{
		const auto min = Point{ DEFAULT_MIN_COORD, DEFAULT_MIN_COORD, DEFAULT_MIN_COORD };
		const auto max = Point{ DEFAULT_MAX_COORD, DEFAULT_MAX_COORD, DEFAULT_MAX_COORD };

		return random_points(n, min, max);
	}

	inline auto random_map(const std::size_t n = 10'000)
	{
		const auto min_coord =
		        get_random(DEFAULT_MIN_COORD - DEFAULT_MAP_VARIANCE, DEFAULT_MIN_COORD + DEFAULT_MAP_VARIANCE);
		const auto max_coord =
		        get_random(DEFAULT_MAX_COORD - DEFAULT_MAP_VARIANCE, DEFAULT_MAX_COORD + DEFAULT_MAP_VARIANCE);

		return random_points(n, { min_coord, min_coord, min_coord }, { max_coord, max_coord, max_coord });
	}

	inline auto naive_box()
	{
		const auto min = Point{ 0.0, 0.0, 0.0 };
		const auto max = Point{ 9.0, 9.0, 9.0 };

		return chs::Box<Point, 3>{ std::make_pair(min, max) };
	}
} // namespace chs::mockup