#pragma once

#include <array>

#include "cheesemap/concepts/Kernel.hpp"
#include "cheesemap/concepts/Point.hpp"

#include "cheesemap/utils/Box.hpp"

namespace chs::kernels
{
	template<typename Point_type, std::size_t Dim = 3>
	        requires chs::concepts::Point<Point_type, Dim>
	class Cube
	{
		Point_type           center_{};
		double               radius_{};
		Box<Point_type, Dim> box_;

		template<std::size_t... Is>
		[[nodiscard]] inline auto is_inside(const auto & p, std::index_sequence<Is...>) const -> bool
		{ return ((box_.min()[Is] <= p[Is] and p[Is] <= box_.max()[Is]) and ...); }

		public:
		static constexpr std::size_t dim = Dim;

		Cube() = delete;
		Cube(const Point_type & center, const double radius) :
		        center_(center), radius_(radius), box_(center_, radius_)
		{}

		[[nodiscard]] inline auto center() const -> const Point_type & { return center_; }
		[[nodiscard]] inline auto radius() const -> double { return radius_; }
		[[nodiscard]] inline auto box() const -> const Box<Point_type, Dim> & { return box_; }

		template<typename QueryPoint>
		        requires chs::concepts::Point<QueryPoint, Dim>
		[[nodiscard]] inline auto is_inside(const QueryPoint & p) const -> bool
		{ return is_inside(p, std::make_index_sequence<Dim>{}); }
	};

	template<typename Point_type>
	Cube(const Point_type &, double) -> Cube<Point_type>;
} // namespace chs::kernels

static_assert(chs::concepts::Kernel<chs::kernels::Cube<std::array<double, 3>>, std::array<double, 3>>);
