#pragma once

#include <array>
#include <unordered_map>
#include <vector>

#include <range/v3/all.hpp>

#include "cheesemap/utils/Box.hpp"
#include "cheesemap/utils/Cell.hpp"

#include "cheesemap/kernels/Sphere.hpp"

#include "cheesemap/concepts/concepts.hpp"

namespace chs
{
	template<typename Point_type, std::size_t Dim = 3>
	class Sparse
	{
		protected:
		using resolution_type   = double;
		using dimensions_array  = std::array<resolution_type, Dim>;
		using dimensions_vector = std::vector<resolution_type>;
		using indices_array     = std::array<std::size_t, Dim>;
		using indices_vector    = std::vector<std::size_t>;
		using cell_type         = Cell<Point_type>;

		static constexpr dimensions_array DEFAULT_RESOLUTIONS = []() {
			std::array<resolution_type, Dim> res{};
			res.fill(1);
			return res;
		}();

		// Dimension of each cell
		dimensions_vector resolutions_{ DEFAULT_RESOLUTIONS.begin(), DEFAULT_RESOLUTIONS.end() };

		// Bounding box of the map
		Box box_;

		// Number of cells of the map on each dimension
		indices_array sizes_{ Dim };

		// Cells of the map
		std::unordered_map<std::size_t, cell_type> cells_;

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

		static auto global2indices(const auto idx, const ranges::range auto & sizes)
		{
			indices_array indices;

			for (const auto i : ranges::views::indices(Dim))
			{
				indices[i] =
				        (idx / ranges::accumulate(sizes | ranges::views::drop(i + 1), std::size_t{ 1 },
				                                  std::multiplies<std::size_t>{})) %
				        sizes[i];
			}

			return indices;
		}

		static auto indices2global(const ranges::range auto & indices, const ranges::range auto & sizes)
		{
			const auto cell_idx = ranges::accumulate(ranges::views::zip(indices, sizes), std::size_t{ 0 },
			                                         [](const auto acc, const auto pair_idx_size) {
				                                         const auto [idx, size] = pair_idx_size;
				                                         return acc * size + idx;
			                                         });
			return cell_idx;
		}

		[[nodiscard]] inline auto indices2global(const ranges::range auto & indices) const
		{
			return indices2global(indices, sizes_);
		}

		[[nodiscard]] inline auto & at(const ranges::range auto & indices)
		{
			return cells_[indices2global(indices, sizes_)];
		}

		[[nodiscard]] inline auto & at(const ranges::range auto & indices) const
		{
			return cells_[indices2global(indices, sizes_)];
		}

		[[nodiscard]] inline auto coord2indices(const Point & p) const
		{
			indices_array idx;
			for (const auto i : ranges::views::indices(Dim))
			{
				const auto rel = p[i] - box_.min()[i];
				idx[i] = static_cast<std::size_t>(std::clamp(rel, 0.0, box_.max()[i] - box_.min()[i]) /
				                                  resolutions_[i]);
			}
			return idx;
		}

		[[nodiscard]] auto submap_dimensions(const ranges::range auto & min,
		                                     const ranges::range auto & max) const
		{
			const auto sizes_view = ranges::views::zip_with(std::minus<>{}, max, min) |
			                        ranges::views::transform([](const auto i) { return i + 1; });
			indices_array sizes;
			ranges::copy(sizes_view, sizes.begin());
			return sizes;
		}

		[[nodiscard]] auto indcs_cells_to_search(const auto & kernel) const
		{
			std::vector<indices_vector> cells_to_search;

			const auto min = coord2indices(kernel.box().min());
			const auto max = coord2indices(kernel.box().max());

			const auto search_dim = submap_dimensions(min, max);

			const auto num_cells =
			        ranges::accumulate(search_dim, std::size_t{ 1 }, std::multiplies<std::size_t>{});

			for (const auto i : ranges::views::indices(num_cells))
			{
				const auto slice_indices  = global2indices(i, search_dim);
				const auto global_indices = ranges::views::zip_with(std::plus<>{}, slice_indices, min);
				cells_to_search.emplace_back(global_indices | ranges::to_vector);
			}

			return cells_to_search;
		}

		public:
		Sparse() = delete;

		template<typename Points_rng>
		Sparse(Points_rng & points, const resolution_type res) : Sparse(points, dimensions_vector(Dim, res))
		{}

		template<typename Points_rng>
		Sparse(Points_rng & points, dimensions_vector res) :
		        resolutions_(std::move(res)), box_(Box::mbb(points))
		{
			for (const auto i : ranges::views::indices(Dim))
			{
				sizes_[i] = static_cast<std::size_t>(
				        std::floor((box_.max()[i] - box_.min()[i]) / resolutions_[i]) + 1);
			}

			for (auto & point : points)
			{
				const auto indices = coord2indices(point);

				const auto global_idx = indices2global(indices, sizes_);

				auto cell_it = cells_.find(global_idx);

				if (cell_it == cells_.end())
				{
					cell_it = cells_.emplace(global_idx, cell_type{ idx2box(indices) }).first;
				}

				auto & cell = cell_it->second;

				cell.add_point(&point);
			}
		}

		template<chs::concepts::Kernel<chs::Point> Kernel_t>
		[[nodiscard]] inline auto query(const Kernel_t & kernel) const
		{
			const auto dummy = []([[maybe_unused]] const auto &) { return true; };
			return query(kernel, dummy);
		}

		template<chs::concepts::Kernel<chs::Point> Kernel_t, chs::concepts::Filter<Point_type> Filter_t>
		[[nodiscard]] inline auto query(const Kernel_t & kernel, Filter_t && filter) const
		{
			std::vector<Point_type *> points;

			const auto min        = coord2indices(kernel.box().min());
			const auto max        = coord2indices(kernel.box().max());
			const auto search_dim = submap_dimensions(min, max);

			const auto num_cells =
			        ranges::accumulate(search_dim, std::size_t{ 1 }, std::multiplies<std::size_t>{});

			for (const auto i : ranges::views::indices(num_cells))
			{
				const auto slice_indices  = global2indices(i, search_dim);
				const auto global_indices = ranges::views::zip_with(std::plus<>{}, slice_indices, min);

				const auto cell_it = cells_.find(indices2global(global_indices, sizes_));

				if (cell_it == cells_.end()) { continue; }

				const auto & cell = cell_it->second;

				for (const auto & point : cell.points())
				{
					if (kernel.is_inside(*point) and filter(*point)) { points.emplace_back(point); }
				}
			}

			return points;
		}

		[[nodiscard]] inline auto knn(const std::integral auto k, const Point_type & p) const
		{
			const auto distance = [&](const Point_type & a, const Point_type & b) {
				return arma::norm(a - b);
			};

			// Store the points and the distance
			using dist_ptr = std::pair<double, Point *>;

			auto cmp_distance = [](const dist_ptr & a, const dist_ptr & b) { return a.first > b.first; };

			std::priority_queue<dist_ptr, std::vector<dist_ptr>, decltype(cmp_distance)> pre_candidates(
			        cmp_distance);

			std::vector<dist_ptr> candidates;

			// Taboo list (to avoid visiting the same cell twice)
			std::set<std::size_t> taboo;

			auto not_visited = [&](const auto & global_idx) { return not taboo.contains(global_idx); };

			// Do an increasing search
			double search_radius = ranges::max(resolutions_);

			const auto num_cells =
			        ranges::accumulate(sizes_, std::size_t{ 1 }, std::multiplies<std::size_t>{});

			while (std::cmp_less(candidates.size(), k) and std::cmp_less(taboo.size(), num_cells))
			{
				// With the new search radius, move pts_and_dist to candidates
				while (not pre_candidates.empty())
				{
					const auto [dist, point] = pre_candidates.top();
					if (dist > search_radius) { break; }
					candidates.emplace_back(dist, point);
					pre_candidates.pop();
				}

				chs::kernels::Sphere<Dim> search(p, search_radius);

				// Get the cells to visit (and haven't been visited yet)
				const auto cells_to_search = indcs_cells_to_search(search);

				auto to_visit = cells_to_search | ranges::views::transform([&](const auto & indices) {
					                return indices2global(indices);
				                }) |
				                ranges::views::filter(not_visited);

				// Visit the cells
				for (const auto cell_idx : to_visit)
				{
					// Mark the cell as visited
					taboo.insert(cell_idx);

					const auto cell_it = cells_.find(cell_idx);

					if (cell_it == cells_.end()) { continue; }

					// Get the cell
					auto & cell = cell_it->second;

					// If the cell is completely inside the search sphere, directly to candidates
					ranges::for_each(cell.points(), [&](const auto & point_ptr) {
						const auto d = distance(p, *point_ptr);
						if (d < search.radius()) { candidates.emplace_back(d, point_ptr); }
						else { pre_candidates.emplace(d, point_ptr); }
					});
				}

				// Update the search radius
				search_radius *= 2.0;
			}

			// Sort the points by distance
			ranges::actions::sort(candidates,
			                      [](const auto & a, const auto & b) { return a.first < b.first; });

			return candidates | ranges::views::take(k) | ranges::to_vector;
		}
	};
} // namespace chs