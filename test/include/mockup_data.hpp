#pragma once

#include <vector>

#include <cheesemap/utils/Point.hpp>

#include "utils.hpp"

namespace chs::mockup
{
	auto naive_points()
	{
		std::vector<chs::Point> points = {
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

	auto random_point(const chs::Point & min, const chs::Point & max)
	{
		chs::Point point;

		for (const auto i : ranges::views::indices(3U))
		{
			point[i] = get_random(min[i], max[i]);
		}

		return point;
	}

	auto random_points(const std::size_t n, const chs::Point & min, const chs::Point & max)
	{
		std::vector<chs::Point> points(n);

		for (auto & point : points)
		{
			point = random_point(min, max);
		}

		return points;
	}
} // namespace chs::mockup