#pragma once

#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/view/indices.hpp>
#include <range/v3/view/transform.hpp>

#include "cheesemap/concepts/Kernel.hpp"

#include "cheesemap/utils/arithmetic.hpp"
#include "cheesemap/utils/Box.hpp"
#include "cheesemap/utils/Point.hpp"

namespace chs::kernels
{
	template<std::size_t Dim = 3>
	class Sphere
	{
		chs::Point center_{};
		double     radius_{};
		double     sq_radius_{};
		chs::Box   box_;

		public:
		Sphere() = delete;
		Sphere(const Point & center, const double radius) :
		        center_(center), radius_(radius), sq_radius_(radius * radius), box_(center_, radius_)
		{}

		[[nodiscard]] inline auto center() const -> const Point & { return center_; }
		[[nodiscard]] inline auto radius() const -> double { return radius_; }
		[[nodiscard]] inline auto box() const -> const Box & { return box_; }

		[[nodiscard]] inline auto is_inside(const Point & p) const -> bool
		{
			return chs::sq_distance<Dim>(center_, p) <= sq_radius_;
		}
	};
} // namespace chs::kernels