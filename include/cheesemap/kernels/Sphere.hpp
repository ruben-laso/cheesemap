#pragma once

#include <array>

#include "cheesemap/concepts/Kernel.hpp"
#include "cheesemap/concepts/Point.hpp"

#include "cheesemap/utils/arithmetic.hpp"
#include "cheesemap/utils/Box.hpp"

namespace chs::kernels
{
	template<typename Point_type, std::size_t Dim = 3>
	        requires chs::concepts::Point<Point_type, Dim>
	class Sphere
	{
		Point_type           center_{};
		double               radius_{};
		double               sq_radius_{};
		Box<Point_type, Dim> box_;

		public:
		static constexpr std::size_t dim = Dim;

		Sphere() = delete;
		Sphere(const Point_type & center, const double radius) :
		        center_(center), radius_(radius), sq_radius_(radius * radius), box_(center_, radius_)
		{}

		[[nodiscard]] inline auto center() const -> const Point_type & { return center_; }
		[[nodiscard]] inline auto radius() const -> double { return radius_; }
		[[nodiscard]] inline auto box() const -> const Box<Point_type, Dim> & { return box_; }

		template<typename QueryPoint>
		        requires chs::concepts::Point<QueryPoint, Dim>
		[[nodiscard]] inline auto is_inside(const QueryPoint & p) const -> bool
		{ return chs::sq_distance<Dim>(center_, p) <= sq_radius_; }
	};

	template<typename Point_type>
	Sphere(const Point_type &, double) -> Sphere<Point_type>;
} // namespace chs::kernels

static_assert(chs::concepts::Kernel<chs::kernels::Sphere<std::array<double, 3>>, std::array<double, 3>>);
