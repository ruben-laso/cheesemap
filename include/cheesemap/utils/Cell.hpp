#pragma once

#include <memory>
#include <vector>

#include <range/v3/all.hpp>

#include "Box.hpp"

namespace chs
{
	template<typename Point_type>
	class Cell
	{
		using Point_ptr_type  = Point_type *;
		using Rng_points_type = std::vector<Point_ptr_type>;

		Box             box_;
		Rng_points_type points_;

		public:
		Cell() = delete;

		explicit Cell(const Box & box) : box_(box) {}

		Cell(const Box & box, ranges::range auto & points) :
		        box_(box),
		        points_(points | ranges::views::transform([](auto & p) { return &p; }) |
		                ranges::to<Rng_points_type>)
		{}

		inline void add_point(Point_ptr_type pnt) { points_.push_back(pnt); }

		[[nodiscard]] inline auto box() const -> const auto & { return box_; }
		// [[nodiscard]] inline auto box() -> auto & { return box_; }

		[[nodiscard]] inline auto points() const -> const auto & { return points_; }
		[[nodiscard]] inline auto points() -> auto & { return points_; }

		[[nodiscard]] inline auto size() const { return points_.size(); }
		[[nodiscard]] inline auto empty() const { return points_.empty(); }
		[[nodiscard]] inline auto begin() const { return points_.begin(); }
		[[nodiscard]] inline auto end() const { return points_.end(); }
	};
} // namespace chs