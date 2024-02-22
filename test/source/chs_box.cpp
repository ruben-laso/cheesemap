#include <cheesemap/utils/Box.hpp>

#include <gtest/gtest.h>

#include "mockup_data.hpp"

TEST(chs, chs_box_build)
{
	const auto points = chs::mockup::naive_points();

	const auto box = chs::Box::mbb(points);

	EXPECT_TRUE(true);
}

TEST(chs, chs_box_build_random)
{
	static constexpr std::size_t num_builds = 100;
	static constexpr std::size_t num_points = 10'000;

	for ([[maybe_unused]] const auto _ : ranges::views::indices(num_builds))
	{
		const auto min_coord = get_random(-10.0, 10.0);
		const auto max_coord = get_random(90.0, 110.0);

		const auto points = chs::mockup::random_points(num_points, { min_coord, min_coord, min_coord },
		                                               { max_coord, max_coord, max_coord });

		EXPECT_NO_THROW(const auto box = chs::Box::mbb(points));
	}
}

TEST(chs, chs_box_max_min)
{
	const auto points = chs::mockup::naive_points();

	const auto box = chs::Box::mbb(points);

	const auto & min = box.min();
	const auto & max = box.max();

	for (const auto i : ranges::views::indices(3U))
	{
		EXPECT_DOUBLE_EQ(min[i], 0.0);
		EXPECT_DOUBLE_EQ(max[i], 9.0);
	}
}

TEST(chs, chs_box_corners)
{
	const auto points = chs::mockup::naive_points();

	const auto box = chs::Box::mbb(points);

	const auto corners = box.corners();

	EXPECT_EQ(corners.size(), 8);

	const auto min = chs::Point{ 0.0, 0.0, 0.0 };
	const auto max = chs::Point{ 9.0, 9.0, 9.0 };

	const auto equal_pts = [](const auto & lhs, const auto & rhs) {
		for (const auto i : ranges::views::indices(3U))
		{
			EXPECT_DOUBLE_EQ(lhs[i], rhs[i]);
		}
	};

	equal_pts(corners[0], min);
	equal_pts(corners[1], max);
	equal_pts(corners[2], chs::Point{ 9.0, 0.0, 0.0 });
	equal_pts(corners[3], chs::Point{ 0.0, 9.0, 9.0 });
	equal_pts(corners[4], chs::Point{ 0.0, 9.0, 0.0 });
	equal_pts(corners[5], chs::Point{ 9.0, 0.0, 9.0 });
	equal_pts(corners[6], chs::Point{ 0.0, 0.0, 9.0 });
	equal_pts(corners[7], chs::Point{ 9.0, 9.0, 0.0 });
}

TEST(chs, chs_box_is_inside)
{
	const auto points = chs::mockup::naive_points();

	const auto box = chs::Box::mbb(points);

	for (const auto & point : points)
	{
		EXPECT_TRUE(box.is_inside(point));
	}

	const auto outside = chs::Point{ 10.0, 10.0, 10.0 };
	EXPECT_FALSE(box.is_inside(outside));
}

TEST(chs, chs_box_is_inside_random)
{
	static constexpr std::size_t num_builds = 100;
	static constexpr std::size_t num_points = 10'000;
	static constexpr std::size_t num_checks = 100;

	for ([[maybe_unused]] const auto _ : ranges::views::indices(num_builds))
	{
		const auto min_coord = get_random(-10.0, 10.0);
		const auto max_coord = get_random(90.0, 110.0);

		const auto points = chs::mockup::random_points(num_points, { min_coord, min_coord, min_coord },
		                                               { max_coord, max_coord, max_coord });

		const auto box = chs::Box::mbb(points);

		for ([[maybe_unused]] const auto i : ranges::views::indices(num_checks))
		{
			const auto point     = chs::mockup::random_point(box.min(), box.max());
			const auto is_inside = box.is_inside(point);
			EXPECT_TRUE(is_inside);
			if (not is_inside)
			{
				std::cout << std::fixed << std::setprecision(2);
				std::cout << "Point: " << point.as_row() << " is not inside the box: {"
				          << box.min().as_row() << " " << box.max().as_row() << "}. Limits: {"
				          << min_coord << " " << max_coord << "}" << '\n';
			}
		}
	}
}

TEST(chs, chs_box_is_outside_random)
{
	static constexpr std::size_t num_builds = 100;
	static constexpr std::size_t num_points = 10'000;
	static constexpr std::size_t num_checks = 100;

	for ([[maybe_unused]] const auto _ : ranges::views::indices(num_builds))
	{
		const auto min_coord = get_random(-10.0, 10.0);
		const auto max_coord = get_random(90.0, 110.0);

		const auto points = chs::mockup::random_points(num_points, { min_coord, min_coord, min_coord },
		                                               { max_coord, max_coord, max_coord });

		const auto box = chs::Box::mbb(points);

		for ([[maybe_unused]] const auto i : ranges::views::indices(num_checks))
		{
			auto outside_point = [&]() -> chs::Point {
				auto       point                = chs::mockup::random_point(box.min(), box.max());
				// Randomly select coordinates to be outside the box (at least one is outside)
				const auto num_coords_to_change = get_random(1U, 3U);
				const auto indices_to_change =
				        ranges::views::indices(3U) | ranges::views::sample(num_coords_to_change);
				// Change the selected coordinates to be outside the box
				for (const auto j : indices_to_change)
				{
					const auto below = get_random(std::numeric_limits<double>::lowest(), min_coord);
					const auto above = get_random(max_coord, std::numeric_limits<double>::max());
					point[j]         = random_bool() ? below : above;
				}
				return point;
			};

			const auto point     = outside_point();
			const auto is_inside = box.is_inside(point);
			EXPECT_FALSE(is_inside);
			if (is_inside)
			{
				std::cout << std::fixed << std::setprecision(2);
				std::cout << "Point: " << point.as_row() << " is inside the box: {"
				          << box.min().as_row() << " " << box.max().as_row() << "}. Limits: {"
				          << min_coord << " " << max_coord << "}" << '\n';
			}
		}
	}
}

auto main() -> int
{
	::testing::InitGoogleTest();
	return RUN_ALL_TESTS();
}