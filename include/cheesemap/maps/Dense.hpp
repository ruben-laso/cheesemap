#pragma once

#include <utility>
#include <vector>

#include <range/v3/all.hpp>

#include "cheesemap/concepts/concepts.hpp"

#include "cheesemap/utils/Box.hpp"
#include "cheesemap/utils/Cell.hpp"

namespace chs
{
	template<typename Point_type, std::size_t Dim = 3>
	class Dense
	{
		using resolution_type = double;
		using array_type      = std::array<resolution_type, Dim>;
		using vector_type     = std::vector<resolution_type>;
		using cell_type       = Cell<Point_type>;

		private:
		static constexpr array_type DEFAULT_RESOLUTIONS = []() {
			std::array<resolution_type, Dim> res{};
			res.fill(1);
			return res;
		}();

		// Dimension of each cell
		vector_type resolutions_{ DEFAULT_RESOLUTIONS };

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
			const auto cell_idx = ranges::accumulate(ranges::views::zip(indices, sizes_), std::size_t{ 0 },
			                                         [](const auto acc, const auto pair_idx_size) {
				                                         const auto [idx, size] = pair_idx_size;
				                                         return acc * size + idx;
			                                         });
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
		Dense() = delete;

		template<typename Points_rng>
		Dense(Points_rng & points, const resolution_type res) : Dense(points, vector_type(Dim, res))
		{}

		template<typename Points_rng>
		Dense(Points_rng & points, const vector_type & res) : resolutions_(res), box_(Box::mbb(points))
		{
			// Number of cells in each dimension
			for (const auto i : ranges::views::indices(Dim))
			{
				sizes_[i] = static_cast<std::size_t>(
				                    std::floor((box_.max()[i] - box_.min()[i]) / resolutions_[i])) +
				            1;
			}

			// Create cells
			const auto num_cells =
			        ranges::accumulate(sizes_, std::size_t{ 1 }, std::multiplies<std::size_t>{});

			for (const auto i : ranges::views::indices(num_cells))
			{
				std::vector<std::size_t> indices(Dim);
				for (const auto j : ranges::views::indices(Dim))
				{
					indices[j] = (i / ranges::accumulate(sizes_ | ranges::views::drop(j + 1),
					                                     std::size_t{ 1 },
					                                     std::multiplies<std::size_t>{})) %
					             sizes_[j];
				}
				cells_.emplace_back(idx2box(indices));
			}

			// Add points to cells
			for (auto & point : points)
			{
				at(coord2idx(point)).add_point(&point);
			}
		}

		template<chs::concepts::Kernel<chs::Point> Kernel_t>
		[[nodiscard]] inline auto query(const Kernel_t & kernel) const
		{
			const auto dummy = []([[maybe_unused]] const auto &) {
				return true;
			};
			return query(kernel, dummy);
		}

		template<chs::concepts::Kernel<chs::Point> Kernel_t, typename Filter_t>
		[[nodiscard]] inline auto query(const Kernel_t & kernel, Filter_t && filter) const
		{
			std::vector<Point *> points;

			std::vector<std::pair<std::size_t, std::size_t>> idxs_min_max(Dim);

			const auto min = coord2idx(kernel.box().min());
			const auto max = coord2idx(kernel.box().max());
			for (const auto i : ranges::views::indices(Dim))
			{
				idxs_min_max[i].first  = std::max(min[i], std::size_t{ 0 });
				idxs_min_max[i].second = std::min(sizes_[i], max[i]);
			}

			std::vector<std::size_t> search_sizes(Dim);
			for (const auto i : ranges::views::indices(Dim))
			{
				search_sizes[i] = idxs_min_max[i].second - idxs_min_max[i].first + 1;
			}

			const auto num_cells =
			        ranges::accumulate(search_sizes, std::size_t{ 1 }, std::multiplies<std::size_t>{});

			for (const auto i : ranges::views::indices(num_cells))
			{
				std::vector<std::size_t> indices(Dim);
				for (const auto j : ranges::views::indices(Dim))
				{
					indices[j] = (i / ranges::accumulate(search_sizes | ranges::views::drop(j + 1),
					                                     std::size_t{ 1 },
					                                     std::multiplies<std::size_t>{})) %
					             search_sizes[j];
					indices[j] += idxs_min_max[j].first;
				}
				const auto & cell = at(indices);
				for (const auto & point : cell.points())
				{
					if (kernel.is_inside(*point) and filter(*point)) { points.emplace_back(point); }
				}
			}

			return points;
		}
	};
} // namespace chs