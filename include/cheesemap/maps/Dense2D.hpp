#pragma once

#include <utility>
#include <vector>

#include <range/v3/all.hpp>

#include "cheesemap/concepts/concepts.hpp"

#include "cheesemap/utils/Box.hpp"
#include "cheesemap/utils/Cell.hpp"

namespace chs
{
	template<typename Point_type>
	class Dense2D
	{
		static constexpr std::size_t Dim = 2;

		using resolution_type = double;
		using array_type      = std::array<resolution_type, Dim>;
		using vector_type     = std::vector<resolution_type>;
		using cell_type       = Cell<Point_type>;

		private:
		static constexpr array_type DEFAULT_RESOLUTIONS = { 1, 1 };

		// Dimension of each cell
		vector_type resolutions_{ DEFAULT_RESOLUTIONS.begin(), DEFAULT_RESOLUTIONS.end() };

		// Bounding box of the map
		Box box_;

		// Number of cells of the map on each dimension
		std::array<std::size_t, Dim> sizes_{ Dim };

		// Cells of the map
		std::vector<cell_type> cells_;

		[[nodiscard]] inline auto idx2box(const ranges::range auto & idx) const
		{
			Point center{};
			Point radii{};

			for (const auto i : ranges::views::indices(Dim))
			{
				center[i] =
				        box_.min()[i] + (static_cast<resolution_type>(idx[i]) + 0.5) * resolutions_[i];
				radii[i] = resolutions_[i] / 2;
			}

			return Box{ center, radii };
		}

		[[nodiscard]] inline auto coord2idx(const Point & p) const
		{
			std::vector<std::size_t> idx(Dim);
			for (const auto i : ranges::views::indices(Dim))
			{
				const auto rel = p[i] - box_.min()[i];
				idx[i] = static_cast<std::size_t>(std::clamp(rel, 0.0, box_.max()[i] - box_.min()[i]) /
				                                  resolutions_[i]);
			}
			return idx;
		}

		[[nodiscard]] inline auto to_cell_idx(const ranges::range auto & indices) const
		{
			const auto cell_idx = indices[0] * sizes_[1] + indices[1];

			return cell_idx;
		}

		[[nodiscard]] inline auto & at(const ranges::range auto & indices)
		{
			return cells_[to_cell_idx(indices)];
		}

		[[nodiscard]] inline auto & at(const ranges::range auto & indices) const
		{
			return cells_[to_cell_idx(indices)];
		}

		public:
		Dense2D() = delete;

		template<typename Points_rng>
		Dense2D(Points_rng & points, const resolution_type res) : Dense2D(points, vector_type(Dim, res))
		{}

		template<typename Points_rng>
		Dense2D(Points_rng & points, vector_type res) : resolutions_(std::move(res)), box_(Box::mbb(points))
		{
			// Number of cells in each dimension
			for (const auto i : ranges::views::indices(Dim))
			{
				sizes_[i] = static_cast<std::size_t>(
				                    std::floor((box_.max()[i] - box_.min()[i]) / resolutions_[i])) +
				            1;
			}

			// Create cells
			cells_.reserve(ranges::accumulate(sizes_, std::size_t{ 1 }, std::multiplies<std::size_t>{}));

			for (const auto i : ranges::views::indices(sizes_[0]))
			{
				for (const auto j : ranges::views::indices(sizes_[1]))
				{
					const auto cell_box = idx2box(std::array{ i, j });
					cells_.emplace_back(cell_box);
				}
			}

			// Add points to cells
			ranges::for_each(points, [this](auto & point) { at(coord2idx(point)).add_point(&point); });
		}

		template<chs::concepts::Kernel<Point_type> Kernel_type>
		[[nodiscard]] inline auto query(const Kernel_type & kernel)
		{
			const auto dummy = []([[maybe_unused]] const auto &) {
				return true;
			};
			return query(kernel, dummy);
		}

		template<chs::concepts::Kernel<Point_type> Kernel_type, chs::concepts::Filter<Point_type> Filter_type>
		[[nodiscard]] inline auto query(const Kernel_type & kernel, Filter_type && filter)
		{
			std::vector<Point *> result;

			// Compute in which cells we have to search
			const auto min = coord2idx(kernel.box().min());
			const auto max = coord2idx(kernel.box().max());

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