#pragma once

#include <array>
#include <vector>

#include <range/v3/all.hpp>

#include "cheesemap/utils/execution.hpp"

#include "cheesemap/utils/Box.hpp"
#include "cheesemap/utils/Cell.hpp"
#include "cheesemap/utils/flags.hpp"
#include "cheesemap/utils/type_traits.hpp"

#include "cheesemap/concepts/concepts.hpp"

#include "cheesemap/maps/SmartSlice.hpp"

namespace chs
{
	template<typename Point_type>
	class Mixed3D
	{
		protected:
		static constexpr std::size_t Dim = 3;

		using resolution_type = double;
		using dimensions_type = chs::type_traits::tuple<resolution_type, Dim>;
		using indices_type    = chs::type_traits::tuple<std::size_t, Dim>;

		static constexpr dimensions_type DEFAULT_RESOLUTIONS = chs::n_tuple<Dim>(resolution_type{ 1 });

		// Dimension of each cell
		dimensions_type resolutions_ = DEFAULT_RESOLUTIONS;

		// Bounding box of the map
		Box box_;

		// Number of cells of the map on each dimension
		indices_type sizes_;

		std::vector<chs::slice::Smart<Point_type>> slices_;

		template<std::size_t... Is>
		[[nodiscard]] inline auto idx2box(const auto & idx, std::index_sequence<Is...>) const
		{
			Point min = box_.min();
			Point max = box_.max();

			(((min[Is] = box_.min()[Is] +
			             static_cast<resolution_type>(std::get<Is>(idx)) * std::get<Is>(resolutions_)),
			  (max[Is] = min[Is] + std::get<Is>(resolutions_))),
			 ...);

			return Box(std::make_pair(min, max));
		}

		[[nodiscard]] inline auto idx2box(const auto & idx) const
		{
			return idx2box(idx, std::make_index_sequence<Dim>{});
		}

		template<std::size_t... Is>
		[[nodiscard]] inline auto coord2indices(const Point & p, std::index_sequence<Is...>) const
		{
			indices_type idx;

			resolution_type diff;

			(((diff = p[Is] - box_.min()[Is]),
			  (diff = std::clamp(diff, 0.0, box_.max()[Is] - box_.min()[Is])),
			  (std::get<Is>(idx) = static_cast<std::size_t>(diff / std::get<Is>(resolutions_)))),
			 ...);

			return idx;
		}

		[[nodiscard]] inline auto coord2indices(const Point & p) const
		{
			return coord2indices(p, std::make_index_sequence<Dim>{});
		}

		template<std::size_t... Is>
		[[nodiscard]] inline auto indices2global(const auto & indices, std::index_sequence<Is...>) const
		{
			std::size_t idx = 0;

			((idx = idx * std::get<Is>(sizes_) + std::get<Is>(indices)), ...);

			return idx;
		}

		[[nodiscard]] inline auto indices2global(const auto & indices) const
		{
			return indices2global(indices, std::make_index_sequence<Dim>{});
		}

		public:
		Mixed3D() = default;

		template<typename Points_rng>
		Mixed3D(Points_rng & points, const resolution_type res, const chs::flags::build::flags_t flags = {}) :
		        Mixed3D(points, dimensions_type(n_tuple<Dim>(res)), flags)
		{}

		template<typename Points_rng>
		Mixed3D(Points_rng & points, dimensions_type res, const chs::flags::build::flags_t flags = {}) :
		        resolutions_(res), box_(Box::mbb(points))
		{
			// Number of cells in each dimension
			[&]<std::size_t... Is>(std::index_sequence<Is...>) {
				((std::get<Is>(sizes_) =
				          static_cast<std::size_t>(std::floor((box_.max()[Is] - box_.min()[Is]) /
				                                              std::get<Is>(resolutions_))) +
				          1),
				 ...);
			}(std::make_index_sequence<Dim>{});

			// Generate the slices
			for (const auto k : ranges::views::indices(std::get<2>(sizes_)))
			{
				const auto box_min = idx2box(indices_type{ 0, 0, k });
				const auto box_max =
				        idx2box(indices_type(std::get<0>(sizes_) - 1, std::get<1>(sizes_) - 1, k + 1));

				Box slice_box{ std::make_pair(box_min.min(), box_max.max()) };
				slices_.emplace_back(slice_box, std::make_tuple(std::get<0>(resolutions_),
				                                                std::get<1>(resolutions_)));
			}

			// Sort points by global idx (should improve locality when querying)
			if (flags & chs::flags::build::REORDER)
			{
				const auto proj = [&](const auto & p) { return indices2global(coord2indices(p)); };
				const auto cmp  = [&](const auto & a, const auto & b) { return proj(a) < proj(b); };
#ifdef __cpp_lib_execution
				if (flags & chs::flags::build::PARALLEL)
				{
					std::sort(std::execution::par_unseq, points.begin(), points.end(), cmp);
				}
				else
#endif
					std::sort(points.begin(), points.end(), cmp);
			}

			// Add the points to the slices
			for (auto & point : points)
			{
				const auto [i, j, k] = coord2indices(point);
				slices_[k].add_point(point);
			}

			if (flags & chs::flags::build::SHRINK_TO_FIT)
			{
				ranges::for_each(slices_, [](auto & slice) { slice.shrink_to_fit(); });
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

			const auto [min_i, min_j, min_k] = coord2indices(kernel.box().min());
			const auto [max_i, max_j, max_k] = coord2indices(kernel.box().max());

			for (const auto k : ranges::views::closed_indices(min_k, max_k))
			{
				auto & slice = slices_[k];

				for (const auto indices :
				     ranges::views::cartesian_product(ranges::views::closed_indices(min_i, max_i),
				                                      ranges::views::closed_indices(min_j, max_j)))
				{
					const auto & cell_opt = slice.at(indices);
					if (not cell_opt.has_value()) { continue; }
					for (auto * point_ptr : cell_opt->get())
					{
						if (kernel.is_inside(*point_ptr) and filter(*point_ptr))
						{
							points.push_back(point_ptr);
						}
					}
				}
			}

			return points;
		}

		[[nodiscard]] inline auto knn(const std::integral auto k_neigh, const Point_type & p) const
		{
			// Store the points and the distance
			chs::sorted_vector<std::pair<double, Point_type *>> candidates(k_neigh);

			// Search radius starts within the cell containing p
			double search_radius = idx2box(coord2indices(p)).distance_to_wall(p, /* inside = */ true);

			auto candidates_within_sq_radius = [&] {
				const auto sq_radius = search_radius * search_radius;
				const auto it        = std::upper_bound(candidates.begin(), candidates.end(), sq_radius,
				                                        [](auto d, auto & pt) { return d < pt.first; });
				return it - candidates.begin();
			};

			// Taboo list (to avoid visiting the same cell twice)
			indices_type taboo_mins;
			indices_type taboo_maxs;

			auto is_tabooed = [&](const indices_type & indices) {
				return chs::within_closed_bounds<Dim>(indices, taboo_mins, taboo_maxs);
			};

			// Do an increasing search
			const double default_radius_increment = chs::min<Dim>(resolutions_);

			// Explore first the neighbors of the cell containing p
			{
				const auto [p_i, p_j, p_k] = coord2indices(p);

				taboo_mins = { p_i, p_j, p_k };
				taboo_maxs = { p_i, p_j, p_k };

				const auto & cell_opt = slices_[p_k].at({ p_i, p_j });

				if (cell_opt.has_value())
				{
					ranges::for_each(cell_opt->get(), [&](const auto & point) {
						candidates.insert({ chs::sq_distance(p, *point), point });
					});
				}
			}

			while (
			        // not enough candidates or last candidate is outside the search radius
			        (std::cmp_less(candidates.size(), k_neigh) or
			         candidates.back().first > (search_radius * search_radius)) and
			        // we have not visited all the cells
			        not chs::all_visited<Dim>(taboo_mins, taboo_maxs, sizes_))
			{
				// Estimate the new required search radius -> k * density -> Saves time ~86% of the queries
				if (not candidates.empty() and search_radius > 0)
				{
					const auto density_based_radius = chs::radius_for_density(
					        candidates_within_sq_radius(), search_radius, k_neigh);
					search_radius = std::min(density_based_radius,
					                         search_radius + default_radius_increment);
				}
				else { search_radius += default_radius_increment; }

				const auto min = coord2indices(p - search_radius);
				const auto max = coord2indices(p + search_radius);

				// If min == taboo_mins and max == taboo_maxs, we have already visited all the cells
				if (chs::all_equal<Dim>(min, taboo_mins) and chs::all_equal<Dim>(max, taboo_maxs))
				{
					continue;
				}

				for (const auto k : ranges::views::closed_indices(std::get<2>(min), std::get<2>(max)))
				{
					const auto & slice = slices_[k];
					for (const auto [i, j] : ranges::views::cartesian_product(
					             ranges::views::closed_indices(std::get<0>(min), std::get<0>(max)),
					             ranges::views::closed_indices(std::get<1>(min), std::get<1>(max))))
					{
						if (is_tabooed({ i, j, k })) { continue; }

						const auto & cell_opt = slice.at({ i, j });

						if (not cell_opt.has_value()) { continue; }

						ranges::for_each(cell_opt->get(), [&](const auto & point) {
							candidates.insert({ chs::sq_distance(p, *point), point });
						});
					}
				}

				taboo_mins = min;
				taboo_maxs = max;
			}

			ranges::for_each(candidates,
			                 [](auto & candidate) { candidate.first = std::sqrt(candidate.first); });

			return candidates;
		}

		[[nodiscard]] inline auto points_per_cell() const
		{
			std::vector<std::size_t> num_points(chs::product<Dim>(sizes_));

			for (const auto & [k, slice] : ranges::views::enumerate(slices_))
			{
				for (const auto [i, j] :
				     ranges::views::cartesian_product(ranges::views::indices(std::get<0>(sizes_)),
				                                      ranges::views::indices(std::get<1>(sizes_))))
				{
					const auto & cell = slice.at({ i, j });
					const auto   idx  = indices2global(std::make_tuple(i, j, k));
					num_points[idx]   = cell.has_value() ? cell->get().size() : 0;
				}
			}

			return num_points;
		}

		[[nodiscard]] inline auto mem_footprint() const
		{
			return ranges::accumulate(slices_, sizeof(*this), [](auto acc, const auto & slice) {
				return acc + slice.mem_footprint();
			});
		}
	};
} // namespace chs