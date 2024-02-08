#pragma once

#include <algorithm>
#include <complex>
#include <numeric>

#include <range/v3/all.hpp>
#include <utility>

#include "Point.hpp"

namespace chs
{
	class Box
	{
		Point center_{};
		Point radii_{};
		Point min_{};
		Point max_{};

		public:
		explicit Box(const Point & center, const double radius) :
		        center_(center),
		        radii_({ radius, radius, radius }),
		        min_(center_ - radii_),
		        max_(center_ + radii_)
		{}

		explicit Box(const Point & center, const Point & radii) :
		        center_(center), radii_(radii), min_(center_ - radii_), max_(center_ + radii_)
		{}

		explicit Box(const std::pair<Point, Point> & min_max) :
		        center_((min_max.first + min_max.second) / 2),
		        radii_(arma::abs(min_max.second - min_max.first) / 2),
		        min_(min_max.first),
		        max_(min_max.second)
		{}

		template<std::size_t Dim = 3>
		[[nodiscard]] inline auto is_inside(const Point & p) const -> bool
		{
			return ranges::all_of(ranges::views::indices(Dim),
			                      [&](const auto i) { return min_[i] <= p[i] && p[i] <= max_[i]; });
		}

		[[nodiscard]] inline auto min() const -> const auto & { return min_; }
		[[nodiscard]] inline auto max() const -> const auto & { return max_; }

		[[nodiscard]] inline auto center() const -> const auto & { return center_; }
		[[nodiscard]] inline auto center() -> auto & { return center_; }
		[[nodiscard]] inline auto radii() const -> const auto & { return radii_; }
		[[nodiscard]] inline auto radii() -> auto & { return radii_; }

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
	};
} // namespace chs