#pragma once

#include <algorithm>
#include <complex>
#include <numeric>

#include <range/v3/all.hpp>
#include <utility>

#include "cheesemap/utils/arithmetic.hpp"
#include "cheesemap/utils/Point.hpp"

namespace chs
{
	class Box
	{
		Point min_{};
		Point max_{};

		public:
		explicit Box(const Point & center, const double radius) :
		        min_(center - Point{ radius, radius, radius }), max_(center + Point{ radius, radius, radius })
		{}

		explicit Box(const Point & center, const Point & radii) : min_(center - radii), max_(center + radii) {}

		explicit Box(const std::pair<Point, Point> & min_max) : min_(min_max.first), max_(min_max.second) {}

		template<std::size_t Dim = 3>
		[[nodiscard]] inline auto is_inside(const Point & p) const -> bool
		{
			return ranges::all_of(ranges::views::indices(Dim),
			                      [&](const auto i) { return min_[i] <= p[i] && p[i] <= max_[i]; });
		}

		[[nodiscard]] inline auto min() const -> const auto & { return min_; }
		[[nodiscard]] inline auto max() const -> const auto & { return max_; }

		[[nodiscard]] inline auto corners() const -> std::vector<Point>
		{
			std::vector<Point> corners;

			corners.push_back(min_);
			corners.push_back(max_);

			// This algorithm works for 3D boxes
			for (const auto i : ranges::views::indices(3U))
			{
				auto min_i = min_;
				auto max_i = max_;

				min_i[i] = max_[i];
				max_i[i] = min_[i];

				corners.push_back(min_i);
				corners.push_back(max_i);
			}

			return corners;
		}

		[[nodiscard]] static auto mbb(const ranges::range auto & points)
		{
			const auto min_max = ranges::accumulate(
			        points,
			        std::pair{ Point{ arma::fill::value(std::numeric_limits<double>::max()) },
			                   Point{ arma::fill::value(std::numeric_limits<double>::lowest()) } },
			        [](const auto & acc, const auto & p) {
				        return std::pair{ arma::min(acc.first, p), arma::max(acc.second, p) };
			        });

			return Box{ min_max };
		}

		[[nodiscard]] auto closest_distance(const Point & p)
		{
			const auto dx = std::min({ chs::distance(p, Point{ min_[0], p[1], p[2] }),
			                           chs::distance(p, Point{ max_[0], p[1], p[2] }) });
			const auto dy = std::min({ chs::distance(p, Point{ p[0], min_[1], p[2] }),
			                           chs::distance(p, Point{ p[0], max_[1], p[2] }) });
			const auto dz = std::min({ chs::distance(p, Point{ p[0], p[1], min_[2] }),
			                           chs::distance(p, Point{ p[0], p[1], max_[2] }) });

			return std::min({ dx, dy, dz });
		}
	};
} // namespace chs