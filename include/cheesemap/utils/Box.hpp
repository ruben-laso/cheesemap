#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>

#include <range/v3/all.hpp>

#include "cheesemap/concepts/Point.hpp"
#include "cheesemap/utils/arithmetic.hpp"

namespace chs
{
	template<typename Point_type, std::size_t Dim = 3>
	        requires chs::concepts::Point<Point_type, Dim>
	class Box
	{
		Point_type min_{};
		Point_type max_{};

		[[nodiscard]] auto closest_distance_within(const Point_type & p) const -> double
		{
			double dist = std::numeric_limits<double>::max();
			for (std::size_t i = 0; i < Dim; ++i)
			{
				const double to_min =
				        std::abs(static_cast<double>(p[i]) - static_cast<double>(min_[i]));
				const double to_max =
				        std::abs(static_cast<double>(p[i]) - static_cast<double>(max_[i]));
				dist = std::min({ dist, to_min, to_max });
			}
			return dist;
		}

		public:
		explicit Box(const Point_type & center, const double radius) : min_(center), max_(center)
		{
			chs::fill<Dim>(min_, [&](const auto i) { return center[i] - radius; });
			chs::fill<Dim>(max_, [&](const auto i) { return center[i] + radius; });
		}

		explicit Box(const Point_type & center, const Point_type & radii) : min_(center), max_(center)
		{
			chs::fill<Dim>(min_, [&](const auto i) { center[i] - radii[i]; });
			chs::fill<Dim>(max_, [&](const auto i) { center[i] + radii[i]; });
		}

		explicit Box(const std::pair<Point_type, Point_type> & min_max) :
		        min_(min_max.first), max_(min_max.second)
		{}

		[[nodiscard]] inline auto is_inside(const Point_type & p) const -> bool
		{
			return ranges::all_of(ranges::views::indices(Dim),
			                      [&](const auto i) { return min_[i] <= p[i] and p[i] <= max_[i]; });
		}

		[[nodiscard]] inline auto min() const -> const auto & { return min_; }
		[[nodiscard]] inline auto max() const -> const auto & { return max_; }

		template<std::size_t Dim_ = Dim>
		[[nodiscard]] inline auto corners() const -> std::array<Point_type, 1 << Dim_>
		{
			std::array<Point_type, 1 << Dim_> corners;

			for (std::size_t i = 0; i < (1 << Dim_); ++i)
			{
				Point_type corner = min_;
				for (std::size_t d = 0; d < Dim_; ++d)
				{
					if ((i >> d) & 1) { corner[d] = max_[d]; }
				}
				corners[i] = corner;
			}

			return corners;
		}

		[[nodiscard]] static auto mbb(const ranges::range auto & points)
		{
			auto lo = points.front();
			auto hi = points.front();

			chs::fill<Dim>(lo, [](const auto) { return std::numeric_limits<double>::max(); });
			chs::fill<Dim>(hi, [](const auto) { return std::numeric_limits<double>::lowest(); });

			for (const auto & p : points)
			{
				for (std::size_t i = 0; i < Dim; ++i)
				{
					lo[i] = std::min<double>(lo[i], p[i]);
					hi[i] = std::max<double>(hi[i], p[i]);
				}
			}

			return Box<Point_type, Dim>{ std::pair{ lo, hi } };
		}

		[[nodiscard]] auto distance(const Point_type & p) const
		{
			auto n = p;
			chs::clamp<Dim>(n, min_, max_);

			return chs::distance<Dim>(p, n);
		}

		[[nodiscard]] auto distance_to_wall(const Point_type & p, const bool inside) const
		{
			if (inside) { return closest_distance_within(p); }
			return distance(p);
		}

		[[nodiscard]] auto distance_to_wall(const Point_type & p) const
		{ return distance_to_wall(p, is_inside(p)); }
	};
} // namespace chs
