#pragma once

#include <algorithm>

#include "cheesemap/concepts/Kernel.hpp"

#include "cheesemap/utils/Box.hpp"
#include "cheesemap/utils/Point.hpp"

namespace chs::kernels
{
	template<std::size_t Dim = 3>
	class Cube
	{
		chs::Point center_{};
		double     radius_{};
		chs::Box   box_;

		public:
		Cube() = delete;
		Cube(const Point & center, const double radius) :
		        center_(center), radius_(radius), box_(center_, radius_)
		{}

		[[nodiscard]] inline auto center() const -> const Point & { return center_; }
		[[nodiscard]] inline auto radius() const -> double { return radius_; }
		[[nodiscard]] inline auto box() const -> const Box & { return box_; }

		[[nodiscard]] inline auto is_inside(const Point & p) const -> bool
		{
			const arma::vec3 diff = p - center_;

			return std::all_of(diff.begin(), diff.begin() + Dim,
			                   [&](const auto x) { return std::abs(x) <= radius_; });
		}
	};
} // namespace chs::kernels