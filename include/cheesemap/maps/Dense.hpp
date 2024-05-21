#pragma once

#include <array>
#include <execution>
#include <utility>
#include <vector>

#include <range/v3/all.hpp>

#include "cheesemap/concepts/concepts.hpp"
#include "cheesemap/kernels/kernels.hpp"

#include "cheesemap/utils/Box.hpp"
#include "cheesemap/utils/Cell.hpp"
#include "cheesemap/utils/sorted_vector.hpp"

#include "cheesemap/utils/arithmetic.hpp"
#include "cheesemap/utils/Cartesian.hpp"
#include <cheesemap/utils/arithmetic.hpp>

namespace chs
{
	template<typename Point_type, std::size_t Dim = 3>
	class Dense
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
		std::vector<cell_type> cells_;

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
			Point min = box_.min();
			Point max = box_.max();

			(((min[Is] = box_.min()[Is] + static_cast<resolution_type>(idx[Is]) * resolutions_[Is]),
			  (max[Is] = min[Is] + resolutions_[Is])),
			 ...);

			return Box(std::make_pair(min, max));
		}

		[[nodiscard]] inline auto idx2box(const ranges::range auto & idx) const
		{
			return idx2box(idx, std::make_index_sequence<Dim>{});
		}

		template<std::size_t... Is>
		[[nodiscard]] inline auto coord2indices(const Point & p, std::index_sequence<Is...>) const
		{
			indices_array idx;

			resolution_type diff;

			(((diff = p[Is] - box_.min()[Is]),
			  (diff = std::clamp(diff, 0.0, box_.max()[Is] - box_.min()[Is])),
			  (idx[Is] = static_cast<std::size_t>(diff / resolutions_[Is]))),
			 ...);

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
			return cells_[indices2global(indices)];
		}

		public:
		Dense() = delete;

		template<ranges::range Points_rng>
		Dense(Points_rng & points, const resolution_type res, const bool reorder = false) :
		        Dense(points, dimensions_vector(Dim, res), reorder)
		{}

		template<ranges::range Points_rng> // requires to be sortable
		Dense(Points_rng & points, dimensions_vector res, const bool reorder = false) :
		        resolutions_(std::move(res)), box_(Box::mbb(points))
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

			cells_.resize(num_cells);

			// Sort points by global idx (should improve locality when querying)

			if (reorder)
			{
				auto proj = [&](const auto & p) { return indices2global(coord2indices(p)); };
				std::sort(std::execution::par_unseq, points.begin(), points.end(),
				          [&](const auto & a, const auto & b) { return proj(a) < proj(b); });
			}

			ranges::for_each(points, [&](auto & point) { at(coord2indices(point)).emplace_back(&point); });
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
			std::vector<Point *> points;

			const auto min = coord2indices(kernel.box().min());
			const auto max = coord2indices(kernel.box().max());

			for (const auto indices : chs::cartesian_as_array<Dim>(min, max))
			{
				const auto & cell = at(indices);
				for (const auto & point : cell)
				{
					if (kernel.is_inside(*point) and filter(*point)) { points.emplace_back(point); }
				}
			}

			return points;
		}

		[[nodiscard]] inline auto knn(const std::integral auto k, const Point_type & p) const
		{
			// Store the points and the distance
			chs::sorted_vector<std::pair<double, Point_type *>> candidates(k);

			// Search radius starts within the cell containing p
			double search_radius = idx2box(coord2indices(p)).distance_to_wall(p, /* inside = */ true);

			auto candidates_within_radius = [&] {
				const auto it = std::upper_bound(candidates.begin(), candidates.end(), search_radius,
				                                 [](auto d, auto & pt) { return d < pt.first; });
				return it - candidates.begin();
			};

			// Taboo list (to avoid visiting the same cell twice)
			indices_array taboo_mins;
			indices_array taboo_maxs;

			auto is_taboo = [&](const auto & indices) {
				return chs::within_closed_bounds<Dim>(indices, taboo_mins, taboo_maxs);
			};

			// Do an increasing search
			const double default_radius_increment = ranges::max(resolutions_);

			// Explore first the neighbors of the cell containing p
			{
				const auto indices = coord2indices(p);
				taboo_mins         = indices;
				taboo_maxs         = indices;
				const auto & cell  = at(indices);
				for (const auto & point : cell)
				{
					candidates.insert({ chs::distance(p, *point), point });
				}
			}

			while (
			        // not enough candidates or last candidate is outside the search radius
			        (std::cmp_less(candidates.size(), k) or candidates.back().first > search_radius) and
			        // we have not visited all the cells
			        not chs::all_visited<Dim>(taboo_mins, taboo_maxs, sizes_))
			{
				// Estimate the new required search radius -> k * density -> Saves time ~86% of the queries
				if (not candidates.empty() and search_radius > 0)
				{
					const auto density_based_radius =
					        chs::radius_for_density(candidates_within_radius(), search_radius, k);
					search_radius = std::min(density_based_radius,
					                         search_radius + default_radius_increment);
				}
				else { search_radius += default_radius_increment; }

				chs::kernels::Sphere<Dim> search(p, search_radius);

				const auto min = coord2indices(search.box().min());
				const auto max = coord2indices(search.box().max());

				// If min == taboo_mins and max == taboo_maxs, we have already visited all the cells
				if (chs::all_equal<Dim>(min, taboo_mins) and chs::all_equal<Dim>(max, taboo_maxs))
				{
					continue;
				}

				for (const auto indices : chs::cartesian_as_array<Dim>(min, max))
				{
					if (is_taboo(indices)) { continue; }

					const auto & cell = at(indices);
					for (const auto & point : cell)
					{
						candidates.insert({ chs::distance(p, *point), point });
					}
				}

				taboo_mins = min;
				taboo_maxs = max;
			}

			return candidates;
		}
	};

} // namespace chs