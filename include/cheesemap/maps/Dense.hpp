#pragma once

#include <array>
#include <execution>
#include <numeric>
#include <queue>
#include <set>
#include <utility>
#include <vector>

#include <range/v3/all.hpp>

#include "cheesemap/concepts/concepts.hpp"
#include "cheesemap/kernels/kernels.hpp"

#include "cheesemap/utils/Box.hpp"
#include "cheesemap/utils/Cell.hpp"

#include "cheesemap/utils/Cartesian.hpp"

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

			// Do an increasing search
			const double default_radius_increment = ranges::max(resolutions_);

			// Explore first the neighbors of the cell containing p
			double search_radius = idx2box(coord2indices(p)).closest_distance(p);
			{
				const auto global_idx = indices2global(coord2indices(p));
				taboo.insert(global_idx);
				const auto & cell = at(coord2indices(p));
				for (const auto & point : cell)
				{
					const auto d = distance(p, *point);
					if (d < search_radius) { candidates.emplace_back(d, point); }
					else { pre_candidates.emplace(d, point); }
				}
			}

			const auto sphere_volume = [](const auto & r) {
				static constexpr auto ratio = 4.0 * std::numbers::pi / 3.0;
				return ratio * std::pow(r, 3);
			};

			while (std::cmp_less(candidates.size(), k) and std::cmp_less(taboo.size(), cells_.size()))
			{
				auto radius_increment = default_radius_increment;
				// Estimate the new required search radius -> k * density
				if (not candidates.empty())
				{
					const auto density =
					        static_cast<double>(candidates.size()) / sphere_volume(search_radius);
					const auto density_based_radius = std::cbrt(static_cast<double>(k) / density);

					if (density_based_radius > search_radius)
					{
						radius_increment = density_based_radius - search_radius;
					}
				}
				search_radius += radius_increment;

				// With the new search radius, move pts_and_dist to candidates
				while (not pre_candidates.empty())
				{
					const auto [dist, point] = pre_candidates.top();
					if (dist > search_radius) { break; }
					candidates.emplace_back(dist, point);
					pre_candidates.pop();
				}

				chs::kernels::Sphere<Dim> search(p, search_radius);

				const auto min = coord2indices(search.box().min());
				const auto max = coord2indices(search.box().max());

				for (const auto indices : chs::cartesian_as_array<Dim>(min, max))
				{
					const auto global_idx = indices2global(indices);
					if (taboo.contains(global_idx)) { continue; }
					else { taboo.insert(global_idx); }

					const auto & cell = cells_[global_idx];
					for (const auto & point : cell)
					{
						const auto d = distance(p, *point);
						if (d < search_radius) { candidates.emplace_back(d, point); }
						else { pre_candidates.emplace(d, point); }
					}
				}
			}

			// Sort the points by distance
			ranges::actions::sort(candidates,
			                      [](const auto & a, const auto & b) { return a.first < b.first; });

			return candidates | ranges::views::take(k) | ranges::to_vector;
		}
	};
} // namespace chs