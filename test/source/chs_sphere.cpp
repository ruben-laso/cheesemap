#include <cheesemap/cheesemap.hpp>

#include <gtest/gtest.h>

#include <range/v3/all.hpp>

#include <armadillo>

#include "mockup_data.hpp"
#include "utils.hpp"

TEST(chs, chs_sphere_build)
{
	chs::Point center{ 0.0, 0.0, 0.0 };

	const auto sphere = chs::kernels::Sphere<3>(center, 1.0);

	for (const auto & i : ranges::views::indices(3U))
	{
		EXPECT_DOUBLE_EQ(center[i], sphere.center()[i]);
	}

	EXPECT_DOUBLE_EQ(1.0, sphere.radius());
}

TEST(chs, chs_sphere_random)
{
	static constexpr std::size_t num_builds = 100;

	for ([[maybe_unused]] const auto _ : ranges::views::indices(num_builds))
	{
		const auto center = chs::mockup::random_point();
		const auto radius = get_random(0.0, 100.0);

		const auto sphere = chs::kernels::Sphere<3>(center, radius);

		for (const auto & i : ranges::views::indices(3U))
		{
			EXPECT_DOUBLE_EQ(center[i], sphere.center()[i]);
		}

		EXPECT_DOUBLE_EQ(radius, sphere.radius());
	}
}

TEST(chs, chs_sphere_inside)
{
	const auto center = chs::Point{ 0.0, 0.0, 0.0 };
	const auto radius = 1.0;

	const auto sphere = chs::kernels::Sphere<3>(center, radius);

	EXPECT_TRUE(sphere.is_inside(chs::Point{ 0.0, 0.0, 0.0 }));
	EXPECT_TRUE(sphere.is_inside(chs::Point{ 0.5, 0.5, 0.5 }));
	EXPECT_TRUE(sphere.is_inside(chs::Point{ 1.0, 0.0, 0.0 }));
	EXPECT_TRUE(sphere.is_inside(chs::Point{ 0.0, 1.0, 0.0 }));
	EXPECT_TRUE(sphere.is_inside(chs::Point{ 0.0, 0.0, 1.0 }));

	EXPECT_FALSE(sphere.is_inside(chs::Point{ 1.0, 1.0, 1.0 }));
	EXPECT_FALSE(sphere.is_inside(chs::Point{ 1.1, 1.1, 1.1 }));
	EXPECT_FALSE(sphere.is_inside(chs::Point{ 2.0, 2.0, 2.0 }));
	EXPECT_FALSE(sphere.is_inside(chs::Point{ 3.0, 3.0, 3.0 }));
}

TEST(chs, chs_sphere_inside_random)
{
	static constexpr std::size_t num_builds = 100;
	static constexpr std::size_t num_points = 100;

	for ([[maybe_unused]] const auto _ : ranges::views::indices(num_builds))
	{
		const auto center = chs::mockup::random_point();
		const auto radius = get_random(0.0, 100.0);

		const auto sphere = chs::kernels::Sphere<3>(center, radius);

		// Generate random points inside the sphere
		for ([[maybe_unused]] const auto & j : ranges::views::indices(num_points))
		{
			const arma::vec3 deviation = arma::normalise<arma::vec3>(
			        { get_random(-1.0, 1.0), get_random(-1.0, 1.0), get_random(-1.0, 1.0) });
			const auto dev_length = get_random(0.0, sphere.radius());

			chs::Point point{ sphere.center() + dev_length * deviation };

			EXPECT_TRUE(sphere.is_inside(point));
			if (not sphere.is_inside(point))
			{
				std::cout << "Point: " << point.as_row()
				          << " is not inside the sphere with center: " << sphere.center().as_row()
				          << " and radius: " << sphere.radius() << '\n';
			}
		}
	}
}

TEST(chs, chs_sphere_outside_random)
{
	static constexpr std::size_t num_builds = 100;
	static constexpr std::size_t num_checks = 100;

	for ([[maybe_unused]] const auto _ : ranges::views::indices(num_builds))
	{
		const auto center = chs::mockup::random_point();
		const auto radius = get_random(0.0, 100.0);

		const auto sphere = chs::kernels::Sphere<3>(center, radius);

		const auto & min = sphere.box().min();
		const auto & max = sphere.box().max();

		for ([[maybe_unused]] const auto & i : ranges::views::indices(num_checks))
		{
			auto outside_point = [&]() -> chs::Point {
				auto       point = chs::mockup::random_point(sphere.box().min(), sphere.box().max());
				// Randomly select coordinates to be outside the box (at least one is outside)
				const auto num_coords_to_change = get_random(1U, 3U);
				const auto indices_to_change =
				        ranges::views::indices(3U) | ranges::views::sample(num_coords_to_change);
				// Change the selected coordinates to be outside the box
				for (const auto j : indices_to_change)
				{
					const auto below = get_random(std::numeric_limits<double>::lowest(), min[j]);
					const auto above = get_random(max[j], std::numeric_limits<double>::max());
					point[j]         = random_bool() ? below : above;
				}
				return point;
			};

			const auto point = outside_point();

			EXPECT_FALSE(sphere.is_inside(point));
		}
	}
}

auto main() -> int
{
	::testing::InitGoogleTest();
	return RUN_ALL_TESTS();
}