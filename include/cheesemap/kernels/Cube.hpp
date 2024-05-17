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

		template<std::size_t... Is>
		[[nodiscard]] inline auto is_inside(const Point & p, std::index_sequence<Is...>) const -> bool
		{
			bool inside = true;

			// This is a fold expression: inside = \forall_i (box.min()[i] <= p[i] <= box.max()[i])
			((inside &= box_.min()[Is] <= p[Is] and p[Is] <= box_.max()[Is]), ...);

			return inside;
		}

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
			return is_inside(p, std::make_index_sequence<Dim>{});
		}
	};
} // namespace chs::kernels