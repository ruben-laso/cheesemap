#pragma once

#include <array>
#include <execution>
#include <queue>
#include <set>
#include <unordered_map>
#include <vector>

#include <range/v3/all.hpp>

#include "cheesemap/utils/Box.hpp"
#include "cheesemap/utils/Cell.hpp"

#include "cheesemap/kernels/Sphere.hpp"

#include "cheesemap/concepts/concepts.hpp"

#include "cheesemap/utils/Cartesian.hpp"
#include "cheesemap/utils/sorted_vector.hpp"

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

		template<std::size_t... Is>
		[[nodiscard]] inline auto indices2global(const ranges::range auto & indices,
		                                         std::index_sequence<Is...>) const
		{
			std::size_t idx = 0;

			((idx = idx * sizes_[Is] + indices[Is]), ...);

			return idx;
		}

		[[nodiscard]] inline auto indices2global(const ranges::range auto & indices) const
		{
			return indices2global(indices, std::make_index_sequence<Dim>{});
		}

		template<std::size_t... Is>
		[[nodiscard]] inline auto idx2box(const ranges::range auto & idx, std::index_sequence<Is...>) const
		{
			Point center{};
			Point radii{};

			((center[Is] =
			          box_.min()[Is] + (static_cast<resolution_type>(idx[Is]) + 0.5) * resolutions_[Is]),
			 ...);
			((radii[Is] = resolutions_[Is] / 2), ...);

			return Box{ center, radii };
		}

		[[nodiscard]] inline auto idx2box(const ranges::range auto & idx) const
		{
			return idx2box(idx, std::make_index_sequence<Dim>{});
		}

		template<std::size_t... Is>
		[[nodiscard]] inline auto coord2indices(const Point & p, std::index_sequence<Is...>) const
		{
			indices_array idx;

			Point diff = p - box_.min();
			((diff[Is] = std::clamp(diff[Is], { 0.0 }, box_.max()[Is] - box_.min()[Is])), ...);

			((idx[Is] = static_cast<std::size_t>(diff[Is] / resolutions_[Is])), ...);

			return idx;
		}

		[[nodiscard]] inline auto coord2indices(const Point & p) const
		{
			return coord2indices(p, std::make_index_sequence<Dim>{});
		}

		[[nodiscard]] inline auto & at(const ranges::range auto & indices)
		{
			return cells_[indices2global(indices)];
		}

		[[nodiscard]] inline auto & at(const ranges::range auto & indices) const
		{
			return cells_.at(indices2global(indices));
		}

		public:
		Sparse() = delete;

		template<typename Points_rng>
		Sparse(Points_rng & points, const resolution_type res, const bool reorder = false) :
		        Sparse(points, dimensions_vector(Dim, res), reorder)
		{}

		template<typename Points_rng>
		Sparse(Points_rng & points, dimensions_vector res, const bool reorder = false) :
		        resolutions_(std::move(res)), box_(Box::mbb(points))
		{
			for (const auto i : ranges::views::indices(Dim))
			{
				sizes_[i] = static_cast<std::size_t>(
				        std::floor((box_.max()[i] - box_.min()[i]) / resolutions_[i]) + 1);
			}

			// Sort points by global idx (should improve locality when querying)
			if (reorder)
			{
				auto proj = [&](const auto & p) { return indices2global(coord2indices(p)); };
				std::sort(std::execution::par_unseq, points.begin(), points.end(),
				          [&](const auto & a, const auto & b) { return proj(a) < proj(b); });
			}

			for (auto & point : points)
			{
				const auto indices = coord2indices(point);

				const auto global_idx = indices2global(indices);

				auto cell_it = cells_.find(global_idx);

				if (cell_it == cells_.end())
				{
					cell_it = cells_.emplace(global_idx, cell_type{}).first;
				}

				auto & cell = cell_it->second;

				cell.emplace_back(&point);
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

			const auto min = coord2indices(kernel.box().min());
			const auto max = coord2indices(kernel.box().max());

			for (const auto indices : chs::cartesian_as_array<Dim>(min, max))
			{
				const auto global_idx = indices2global(indices);
				const auto cell_it    = cells_.find(global_idx);

				if (cell_it == cells_.end()) { continue; }

				const auto & cell = cell_it->second;

				for (const auto & point : cell)
				{
					if (kernel.is_inside(*point) and filter(*point)) { points.emplace_back(point); }
				}
			}

			return points;
		}

		[[nodiscard]] inline auto knn(const std::integral auto k, const Point_type & p) const
		{
			static constexpr auto max_growth_factor = 1.0; // n times the default increment

			const auto distance = [&](const Point_type & a, const Point_type & b) {
				double dist = 0.0;
				dist += (a[0] - b[0]) * (a[0] - b[0]);
				dist += (a[1] - b[1]) * (a[1] - b[1]);
				dist += (a[2] - b[2]) * (a[2] - b[2]);
				return std::sqrt(dist);
			};

			// Store the points and the distance
			chs::sorted_vector<std::pair<double, Point_type *>> candidates(k);

			// Search radius starts within the cell containing p
			double search_radius = idx2box(coord2indices(p)).closest_distance(p);

			auto candidates_within_radius = [&] {
				const auto it = std::upper_bound(candidates.begin(), candidates.end(), search_radius,
				                                 [](auto d, auto & pt) { return d < pt.first; });
				return it - candidates.begin();
			};

			// Taboo list (to avoid visiting the same cell twice)
			indices_array taboo_mins;
			indices_array taboo_maxs;

			auto is_tabooed = [&](const auto & indices) {
				auto within_bounds = [](const auto & tuple_idx_min_max) {
					const auto & [idx, min, max] = tuple_idx_min_max;
					return std::cmp_greater_equal(idx, min) and std::cmp_less_equal(idx, max);
				};
				return ranges::all_of(ranges::views::zip(indices, taboo_mins, taboo_maxs),
				                      within_bounds);
			};

			// Do an increasing search
			const double default_radius_increment = ranges::max(resolutions_);

			// Explore first the neighbors of the cell containing p
			{
				const auto indices = coord2indices(p);
				taboo_mins         = indices;
				taboo_maxs         = indices;

				const auto global_idx = indices2global(indices);
				const auto cell_it    = cells_.find(global_idx);

				if (cell_it != cells_.end())
				{
					const auto & cell = cell_it->second;
					for (const auto & point : cell)
					{
						const auto d = distance(p, *point);
						candidates.insert({ d, point });
					}
				}
			}

			const auto sphere_volume = [](const auto & r) {
				static constexpr auto ratio = 4.0 * std::numbers::pi / 3.0;
				return ratio * std::pow(r, 3);
			};

			while (
			        // not enough candidates or last candidate is outside the search radius
			        (std::cmp_less(candidates.size(), k) or candidates.back().first > search_radius) and
			        // we have not visited all the cells
			        ranges::any_of(ranges::views::zip(taboo_mins, taboo_maxs, sizes_), [](const auto & t) {
				        const auto & [min, max, size] = t;
				        return std::cmp_less(max - min, size - 1);
			        }))
			{
				auto radius_increment = default_radius_increment;
				// Estimate the new required search radius -> k * density -> Saves time ~86% of the queries
				if (not candidates.empty() and search_radius > 0)
				{
					const auto density = static_cast<double>(candidates_within_radius()) /
					                     sphere_volume(search_radius);
					const auto density_based_radius = std::cbrt(static_cast<double>(k) / density);

					radius_increment = std::min(density_based_radius - search_radius,
					                            max_growth_factor * default_radius_increment);
				}
				search_radius += radius_increment;

				chs::kernels::Sphere<Dim> search(p, search_radius);

				const auto min = coord2indices(search.box().min());
				const auto max = coord2indices(search.box().max());

				for (const auto indices : chs::cartesian_as_array<Dim>(min, max))
				{
					if (is_tabooed(indices)) { continue; }

					const auto global_idx = indices2global(indices);

					const auto cell_it = cells_.find(global_idx);

					if (cell_it == cells_.end()) { continue; }

					const auto & cell = cell_it->second;

					for (const auto & point : cell)
					{
						const auto d = distance(p, *point);
						candidates.insert({ d, point });
					}
				}

				taboo_mins = min;
				taboo_maxs = max;
			}

			return candidates;
		}
	};
} // namespace chs