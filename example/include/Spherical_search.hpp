#pragma once

#include <range/v3/all.hpp>

#include <cellmap/concepts/Kernel.hpp>

namespace chs
{
	template<std::size_t Dim = 3>
	class Spherical_search
	{
		chs::Point center_{};
		chs::Box   box_;
		double     radius_{ 1.0 };

		public:
		Spherical_search() = delete;
		Spherical_search(const Point & center, const double radius) :
		        center_(center), box_(center_, radius), radius_(radius)
		{}

		[[nodiscard]] inline auto center() const -> const Point & { return center_; }
		[[nodiscard]] inline auto radius() const -> double { return radius_; }
		[[nodiscard]] inline auto box() const -> const Box & { return box_; }

		[[nodiscard]] inline auto is_inside(const Point & p) const -> bool
		{
			const arma::vec3 diff = p - center_;

			// Norm in Dim dimensions
			auto norm_view = ranges::views::indices(Dim) |
			                 ranges::views::transform([&](const auto i) { return diff[i]; }) |
			                 ranges::views::transform([](const auto x) { return x * x; });

			const auto norm = ranges::accumulate(norm_view, 0.0);

			return norm <= radius_ * radius_;
		}
	};
} // namespace chs

static_assert(chs::concepts::Kernel<chs::Spherical_search<3>, chs::Point>);