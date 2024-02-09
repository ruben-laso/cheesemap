#include <cheesemap/cheesemap.hpp>

#include <gtest/gtest.h>

#include "mockup_data.hpp"

TEST(chs, chs_cell_build)
{
	const auto box = chs::mockup::naive_box();

	EXPECT_NO_THROW(const auto cell = chs::Cell<chs::Point>(box));
}

TEST(chs, chs_cell_build_points)
{
	auto       points = chs::mockup::naive_points();
	const auto box    = chs::Box::mbb(points);

	EXPECT_NO_THROW(const auto cell = chs::Cell<chs::Point>(box, points));
}

TEST(chs, chs_cell_build_random)
{
	static constexpr std::size_t num_builds = 100;
	static constexpr std::size_t num_points = 10'000;

	for ([[maybe_unused]] const auto _ : ranges::views::indices(num_builds))
	{
		const auto min_coord = get_random(-10.0, 10.0);
		const auto max_coord = get_random(90.0, 110.0);

		auto points = chs::mockup::random_points(num_points, { min_coord, min_coord, min_coord },
		                                         { max_coord, max_coord, max_coord });

		const auto box = chs::Box::mbb(points);

		EXPECT_NO_THROW(const auto cell = chs::Cell<chs::Point>(box, points));
	}
}

TEST(chs, chs_cell_points_random)
{
	static constexpr std::size_t num_builds = 100;
	static constexpr std::size_t num_points = 10'000;

	for ([[maybe_unused]] const auto _ : ranges::views::indices(num_builds))
	{
		const auto min_coord = get_random(-10.0, 10.0);
		const auto max_coord = get_random(90.0, 110.0);

		auto points = chs::mockup::random_points(num_points, { min_coord, min_coord, min_coord },
		                                         { max_coord, max_coord, max_coord });

		const auto box  = chs::Box::mbb(points);
		const auto cell = chs::Cell<chs::Point>(box, points);

		EXPECT_TRUE(
		        are_the_same(points | ranges::views::transform([](auto & p) { return &p; }), cell.points()));
	}
}

int main(int argc, char ** argv)
{
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}