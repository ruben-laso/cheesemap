#include <array>
#include <vector>

#include <armadillo>
#include <Eigen/Core>
#include <gtest/gtest.h>

#include "cheesemap/concepts/Point.hpp"
#include "cheesemap/kernels/Cube.hpp"
#include "cheesemap/kernels/Sphere.hpp"
#include "cheesemap/maps/Dense.hpp"
#include "cheesemap/maps/Mixed.hpp"
#include "cheesemap/maps/Sparse.hpp"
#include "cheesemap/utils/Box.hpp"

namespace
{
	struct NotAPoint
	{};

	struct KernelWithoutBox
	{
		auto is_inside(const std::array<double, 3> &) const -> bool { return true; }
	};

	static_assert(not chs::concepts::Point<NotAPoint, 3>);
	static_assert(chs::concepts::Point<arma::vec3, 3>);
	static_assert(chs::concepts::Point<Eigen::Vector3d, 3>);
	static_assert(chs::concepts::Point<std::array<double, 3>, 3>);
	static_assert(chs::concepts::Point<std::vector<double>, 3>);
	static_assert(chs::concepts::Kernel<chs::kernels::Cube<std::array<double, 3>, 3>, std::array<double, 3>>);
	static_assert(chs::concepts::Kernel<chs::kernels::Cube<std::array<double, 3>, 3>, std::vector<double>>);
	static_assert(not chs::concepts::Kernel<KernelWithoutBox, std::array<double, 3>>);

	template<typename Point_type, typename... Coordinates>
	auto make_point(Coordinates... coordinates) -> Point_type
	{
		const std::array values{ static_cast<double>(coordinates)... };
		Point_type       point{};
		if constexpr (requires {
			              point.size();
			              point.resize(values.size());
		              })
		{
			if (point.size() == 0) { point.resize(values.size()); }
		}
		for (std::size_t i = 0; i < values.size(); ++i)
		{
			point[i] = values[i];
		}
		return point;
	}

	template<typename Point_type>
	auto make_points() -> std::vector<Point_type>
	{
		return { make_point<Point_type>(0.0, 0.0, 0.0), make_point<Point_type>(1.0, 1.0, 1.0),
			 make_point<Point_type>(3.0, 3.0, 3.0) };
	}

	template<typename Map_type, typename Point_type>
	void expect_map_query(Map_type & map)
	{
		const auto sphere  = chs::kernels::Sphere<Point_type, 3>(make_point<Point_type>(1.0, 1.0, 1.0), 0.1);
		const auto results = map.query(sphere);

		ASSERT_EQ(results.size(), 1);
		for (std::size_t i = 0; i < 3; ++i)
		{
			EXPECT_DOUBLE_EQ((*results.front())[i], 1.0);
		}
	}

	template<typename Point_type>
	class ThreeDimensionalPointTypes : public ::testing::Test
	{};

	using Point_types = ::testing::Types<arma::vec3, Eigen::Vector3d, std::array<double, 3>>;
	TYPED_TEST_SUITE(ThreeDimensionalPointTypes, Point_types);

	TYPED_TEST(ThreeDimensionalPointTypes, BoxAndSphereAcceptPointType)
	{
		const auto points = make_points<TypeParam>();
		const auto box    = chs::Box<TypeParam, 3>::mbb(points);

		for (std::size_t i = 0; i < 3; ++i)
		{
			EXPECT_DOUBLE_EQ(box.min()[i], 0.0);
			EXPECT_DOUBLE_EQ(box.max()[i], 3.0);
		}

		const auto sphere = chs::kernels::Sphere<TypeParam, 3>(make_point<TypeParam>(1.0, 1.0, 1.0), 0.1);
		EXPECT_TRUE(sphere.is_inside(make_point<TypeParam>(1.0, 1.0, 1.0)));
		EXPECT_FALSE(sphere.is_inside(make_point<TypeParam>(2.0, 2.0, 2.0)));

		const auto cube = chs::kernels::Cube<TypeParam, 3>(make_point<TypeParam>(1.0, 1.0, 1.0), 0.1);
		EXPECT_TRUE(cube.is_inside(make_point<TypeParam>(1.0, 1.0, 1.0)));
		EXPECT_FALSE(cube.is_inside(make_point<TypeParam>(2.0, 2.0, 2.0)));
	}

	TYPED_TEST(ThreeDimensionalPointTypes, DenseSparseAndMixedMapsQuery)
	{
		auto points = make_points<TypeParam>();

		chs::Dense<TypeParam, 3> dense(points, 1.0);
		expect_map_query<decltype(dense), TypeParam>(dense);

		chs::Sparse<TypeParam, 3> sparse(points, 1.0);
		expect_map_query<decltype(sparse), TypeParam>(sparse);

		chs::Mixed<TypeParam, 3> mixed(points, 1.0);
		expect_map_query<decltype(mixed), TypeParam>(mixed);
	}

	TYPED_TEST(ThreeDimensionalPointTypes, PointDimensionKernelQueriesDifferentGridDimension)
	{
		const auto make_points = [] {
			return std::vector<TypeParam>{ make_point<TypeParam>(1.0, 1.0, 0.0),
				                       make_point<TypeParam>(1.0, 1.0, 1.0),
				                       make_point<TypeParam>(1.0, 1.0, 2.0),
				                       make_point<TypeParam>(1.2, 1.2, 2.0) };
		};
		const auto cube = chs::kernels::Cube<TypeParam, 3>(make_point<TypeParam>(1.0, 1.0, 1.0), 0.1);

		auto                     points = make_points();
		chs::Mixed<TypeParam, 2> mixed(points, 1.0);
		const auto               results = mixed.query(cube);
		ASSERT_EQ(results.size(), 1);
		EXPECT_DOUBLE_EQ((*results.front())[0], 1.0);
		EXPECT_DOUBLE_EQ((*results.front())[1], 1.0);
		EXPECT_DOUBLE_EQ((*results.front())[2], 1.0);
	}

	TEST(HeterogeneousPointTypes, ArrayKernelQueriesVectorPointMap)
	{
		using MapPoint    = std::vector<double>;
		using KernelPoint = std::array<double, 3>;

		auto points = std::vector<MapPoint>{
			MapPoint{ 1.0, 1.0, 1.0 },
			MapPoint{ 1.0, 1.0, 0.0 },
			MapPoint{ 1.0, 1.0, 2.0 },
			MapPoint{ 2.0, 2.0, 2.0 },
		};
		chs::Dense<MapPoint, 3> map(points, 1.0);
		const auto cube = chs::kernels::Cube<KernelPoint, 3>(make_point<KernelPoint>(1.0, 1.0, 1.0), 0.1);

		const auto results = map.query(cube);

		ASSERT_EQ(results.size(), 1);
		EXPECT_DOUBLE_EQ((*results.front())[0], 1.0);
		EXPECT_DOUBLE_EQ((*results.front())[1], 1.0);
		EXPECT_DOUBLE_EQ((*results.front())[2], 1.0);

		const auto sphere = chs::kernels::Sphere<KernelPoint, 3>(make_point<KernelPoint>(1.0, 1.0, 1.0), 0.1);
		const auto sphere_results = map.query(sphere);
		ASSERT_EQ(sphere_results.size(), 1);
	}

	template<typename Point_type>
	class TwoDimensionalPointTypes : public ::testing::Test
	{};

	using TwoDimensionalPoint_types = ::testing::Types<arma::vec3, Eigen::Vector3d, std::array<double, 3>>;
	TYPED_TEST_SUITE(TwoDimensionalPointTypes, TwoDimensionalPoint_types);

	TYPED_TEST(TwoDimensionalPointTypes, Mixed2DMapQueriesPointType)
	{
		auto                     points = std::vector<TypeParam>{ make_point<TypeParam>(0.0, 0.0, 0.0),
			                                                  make_point<TypeParam>(1.0, 1.0, 0.0),
			                                                  make_point<TypeParam>(2.25, 2.25, 0.0) };
		chs::Mixed<TypeParam, 2> mixed(points, 1.0);

		const auto sphere  = chs::kernels::Sphere<TypeParam, 2>(make_point<TypeParam>(1.0, 1.0, 0.0), 0.1);
		const auto results = mixed.query(sphere);

		ASSERT_EQ(results.size(), 1);
		EXPECT_DOUBLE_EQ((*results.front())[0], 1.0);
		EXPECT_DOUBLE_EQ((*results.front())[1], 1.0);
	}
} // namespace

auto main() -> int
{
	::testing::InitGoogleTest();
	return RUN_ALL_TESTS();
}
