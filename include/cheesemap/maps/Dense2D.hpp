#pragma once

#include <utility>
#include <vector>

#include <range/v3/all.hpp>

#include "cheesemap/maps/Dense.hpp"

#include "cheesemap/concepts/concepts.hpp"

#include "cheesemap/utils/Box.hpp"
#include "cheesemap/utils/Cell.hpp"

namespace chs
{
	template<typename Point_type>
	class Dense2D : public chs::Dense<Point_type, 2>
	{
		protected:
		using super_type = chs::Dense<Point_type, 2>;

		[[nodiscard]] inline auto & at(const ranges::range auto & indices)
		{
			return this->cells_[indices[0] * this->sizes_[1] + indices[1]];
		}

		[[nodiscard]] inline auto & at(const ranges::range auto & indices) const
		{
			return this->cells_[indices[0] * this->sizes_[1] + indices[1]];
		}

		public:
		Dense2D() = delete;

		template<typename... Args_types>
		Dense2D(Args_types &&... args) : super_type(std::forward<Args_types>(args)...)
		{}

		template<chs::concepts::Kernel<Point_type> Kernel_type>
		[[nodiscard]] inline auto query(const Kernel_type & kernel) const
		{
			const auto dummy = []([[maybe_unused]] const auto &) { return true; };
			return query(kernel, dummy);
		}

		template<chs::concepts::Kernel<Point_type> Kernel_type, chs::concepts::Filter<Point_type> Filter_type>
		[[nodiscard]] inline auto query(const Kernel_type & kernel, Filter_type && filter) const
		{
			std::vector<Point *> result;

			// Compute in which cells we have to search
			const auto min = this->coord2indices(kernel.box().min());
			const auto max = this->coord2indices(kernel.box().max());

			for (const auto i : ranges::views::closed_indices(min[0], max[0]))
			{
				for (const auto j : ranges::views::closed_indices(min[1], max[1]))
				{
					auto & cell = at(std::array{ i, j });
					for (auto * point_ptr : cell.points())
					{
						if (kernel.is_inside(*point_ptr) and filter(*point_ptr))
						{
							result.push_back(point_ptr);
						}
					}
				}
			}

			return result;
		}
	};
} // namespace chs