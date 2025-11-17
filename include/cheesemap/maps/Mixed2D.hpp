#pragma once

#include <array>
#include <vector>

#include <range/v3/all.hpp>

#include "cheesemap/utils/execution.hpp"

#include "cheesemap/concepts/concepts.hpp"

#include "cheesemap/maps/SmartSlice.hpp"

#include "cheesemap/utils/Cartesian.hpp"
#include "cheesemap/utils/flags.hpp"
#include "cheesemap/utils/sorted_vector.hpp"
#include "cheesemap/utils/type_traits.hpp"

namespace chs
{
	template<typename Point_type, template<typename, typename...> class HashMap = std::unordered_map>
	class Mixed2D
	{
		protected:
		static constexpr std::size_t Dim = 2;

		using resolution_type = double;
		using dimensions_type = chs::type_traits::tuple<resolution_type, Dim>;
		using indices_type    = chs::type_traits::tuple<std::size_t, Dim>;

		static constexpr dimensions_type DEFAULT_RESOLUTIONS = chs::n_tuple<Dim>(resolution_type{ 1 });

		chs::slice::Smart<Point_type, HashMap> slice_;

		public:
		Mixed2D() = default;

		template<typename Points_rng>
		Mixed2D(Points_rng & points, const resolution_type res, const chs::flags::build::flags_t flags = {}) :
		        Mixed2D(points, dimensions_type(chs::n_tuple<Dim>(res)), flags)
		{}

		template<typename Points_rng>
		Mixed2D(Points_rng & points, const dimensions_type & res, const chs::flags::build::flags_t flags = {}) :
		        slice_(Box::mbb(points), res)
		{
			// Sort points by global idx (should improve locality when querying)
			if (flags & chs::flags::build::REORDER)
			{
				const auto proj = [&](const auto & p) {
					return slice_.indices2global(slice_.coord2indices(p));
				};
				const auto cmp = [&](const auto & a, const auto & b) { return proj(a) < proj(b); };
#ifdef __cpp_lib_execution
				if (flags & chs::flags::build::PARALLEL)
				{
					std::sort(std::execution::par_unseq, points.begin(), points.end(), cmp);
				}
				else
#endif
					std::sort(points.begin(), points.end(), cmp);
			}

			for (auto & point : points)
			{
				slice_.add_point(point);
			}

			if (flags & chs::flags::build::SHRINK_TO_FIT) { slice_.shrink_to_fit(); }
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

			const auto min = slice_.coord2indices(kernel.box().min());
			const auto max = slice_.coord2indices(kernel.box().max());

			for (const auto indices : chs::cartesian<Dim>(min, max))
			{
				const auto & cell_opt = slice_.at(indices);
				if (not cell_opt.has_value()) { continue; }
				for (auto * point_ptr : cell_opt->get())
				{
					if (kernel.is_inside(*point_ptr) and filter(*point_ptr))
					{
						points.push_back(point_ptr);
					}
				}
			}

			return points;
		}

		[[nodiscard]] inline auto knn(const std::integral auto k, const Point_type & p) const
		{
			// Store the points and the distance
			chs::sorted_vector<std::pair<double, Point_type *>> candidates(k);

			// Search radius starts within the cell containing p
			double search_radius =
			        slice_.idx2box(slice_.coord2indices(p)).distance_to_wall(p, /* inside = */ true);

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
			const double default_radius_increment = chs::min<Dim>(slice_.resolutions());

			// Explore first the neighbors of the cell containing p
			{
				const auto indices = slice_.coord2indices(p);
				taboo_mins         = indices;
				taboo_maxs         = indices;

				const auto & cell_opt = slice_.at(indices);

				if (cell_opt.has_value())
				{
					ranges::for_each(cell_opt->get(), [&](const auto & point) {
						candidates.insert({ chs::sq_distance(p, *point), point });
					});
				}
			}

			while (
			        // not enough candidates or last candidate is outside the search radius
			        (std::cmp_less(candidates.size(), k) or
			         candidates.back().first > (search_radius * search_radius)) and
			        // we have not visited all the cells
			        not chs::all_visited<Dim>(taboo_mins, taboo_maxs, slice_.sizes()))
			{
				// Estimate the new required search radius -> k * density -> Saves time ~86% of the queries
				if (not candidates.empty() and search_radius > 0)
				{
					const auto density_based_radius = chs::radius_for_density(
					        candidates_within_sq_radius(), search_radius, k);
					search_radius = std::min(density_based_radius,
					                         search_radius + default_radius_increment);
				}
				else
				{
					search_radius += default_radius_increment;
				}

				const auto min = slice_.coord2indices(p - search_radius);
				const auto max = slice_.coord2indices(p + search_radius);

				// If min == taboo_mins and max == taboo_maxs, we have already visited all the cells
				if (chs::all_equal<Dim>(min, taboo_mins) and chs::all_equal<Dim>(max, taboo_maxs))
				{
					continue;
				}

				for (const auto indices : chs::cartesian<Dim>(min, max))
				{
					if (is_tabooed(indices)) { continue; }

					const auto & cell_opt = slice_.at(indices);

					if (not cell_opt.has_value()) { continue; }

					ranges::for_each(cell_opt->get(), [&](const auto & point) {
						candidates.insert({ chs::sq_distance(p, *point), point });
					});
				}

				taboo_mins = min;
				taboo_maxs = max;
			}

			ranges::for_each(candidates,
			                 [](auto & candidate) { candidate.first = std::sqrt(candidate.first); });

			return candidates;
		}

		[[nodiscard]] inline auto cells_stored() const
		{
			std::vector<bool> cells_stored(chs::product<Dim>(slice_.sizes()), false);
			for (const auto indices :
			     ranges::views::cartesian_product(ranges::views::indices(std::get<0>(slice_.sizes())),
			                                      ranges::views::indices(std::get<1>(slice_.sizes()))))

			{
				const auto & cell = slice_.at(indices);
				const auto   idx  = slice_.indices2global(indices);
				cells_stored[idx] = cell.has_value();
			}
			return cells_stored;
		}

		[[nodiscard]] inline auto points_per_cell() const
		{
			const auto sizes = slice_.sizes();

			std::vector<std::size_t> num_points(std::get<0>(sizes) * std::get<1>(sizes));

			for (const auto indices :
			     ranges::views::cartesian_product(ranges::views::indices(std::get<0>(sizes)),
			                                      ranges::views::indices(std::get<1>(sizes))))

			{
				const auto & cell = slice_.at(indices);
				const auto   idx  = slice_.indices2global(indices);
				num_points[idx]   = cell.has_value() ? cell->get().size() : 0;
			}

			return num_points;
		}
	};
} // namespace chs