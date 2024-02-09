#pragma once

#include <cheesemap/concepts/Kernel.hpp>

#include <cheesemap/utils/Box.hpp>
#include <cheesemap/utils/Point.hpp>

namespace chs::kernels
{
	template<std::size_t Dim = 3>
	class Sphere
	{
		chs::Point center_{};
		double     radius_{};
		chs::Box   box_;

		public:
		Sphere() = delete;
		Sphere(const Point & center, const double radius) :
		        center_(center), radius_(radius), box_(center_, radius_)
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
} // namespace chs::kernels