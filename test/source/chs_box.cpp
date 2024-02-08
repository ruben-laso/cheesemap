#include <cheesemap/utils/Box.hpp>

#include <gtest/gtest.h>

auto mockup_naive_points()
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

TEST(chs, chs_box_build)
{
	const auto points = mockup_naive_points();

	const auto box = chs::Box::mbb(points);

	EXPECT_TRUE(true);
}

TEST(chs, chs_box_max_min)
{
	const auto points = mockup_naive_points();

	const auto box = chs::Box::mbb(points);

	const auto & min = box.min();
	const auto & max = box.max();

	for (const auto i : ranges::views::indices(3U))
	{
		EXPECT_DOUBLE_EQ(min[i], 0.0);
		EXPECT_DOUBLE_EQ(max[i], 9.0);
	}
}

TEST(chs, chs_box_center)
{
	const auto points = mockup_naive_points();

	const auto box = chs::Box::mbb(points);

	const auto & center = box.center();

	const auto midpoint = chs::Point{ 4.5, 4.5, 4.5 };

	for (const auto i : ranges::views::indices(3U))
	{
		EXPECT_DOUBLE_EQ(center[i], midpoint[i]);
	}
}

TEST(chs, chs_box_radii)
{
	const auto points = mockup_naive_points();

	const auto box = chs::Box::mbb(points);

	const auto & radii = box.radii();

	const auto half = chs::Point{ 4.5, 4.5, 4.5 };

	for (const auto i : ranges::views::indices(3U))
	{
		EXPECT_DOUBLE_EQ(radii[i], half[i]);
	}
}

TEST(chs, chs_box_corners)
{
	const auto points = mockup_naive_points();

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
	const auto points = mockup_naive_points();

	const auto box = chs::Box::mbb(points);

	for (const auto & point : points)
	{
		EXPECT_TRUE(box.is_inside(point));
	}

	const auto outside = chs::Point{ 10.0, 10.0, 10.0 };
	EXPECT_FALSE(box.is_inside(outside));
}


auto main() -> int
{
	::testing::InitGoogleTest();
	return RUN_ALL_TESTS();
}