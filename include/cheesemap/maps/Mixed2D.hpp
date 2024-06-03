#pragma once

#include <execution>
#include <vector>

#include <range/v3/all.hpp>

#include "cheesemap/concepts/concepts.hpp"

#include "cheesemap/kernels/kernels.hpp"

#include "cheesemap/maps/SmartSlice.hpp"

#include "cheesemap/utils/flags.hpp"
#include "cheesemap/utils/sorted_vector.hpp"

namespace chs
{
	template<typename Point_type>
	class Mixed2D
	{
		protected:
		static constexpr std::size_t Dim = 2;

		using resolution_type   = double;
		using dimensions_vector = std::vector<resolution_type>;
		using indices_array     = std::array<std::size_t, Dim>;

		chs::slice::Smart<Point_type> slice_;

		public:
		Mixed2D() = default;

		template<typename Points_rng>
		Mixed2D(Points_rng & points, const resolution_type res,
		        const chs::flags::build::flags_t flags = {}) :
		        Mixed2D(points, dimensions_vector(Dim, res), flags)
		{}

		template<typename Points_rng>
		Mixed2D(Points_rng & points, const dimensions_vector & res,
		        const chs::flags::build::flags_t flags = {}) :
		        slice_(Box::mbb(points), res)
		{
			// Sort points by global idx (should improve locality when querying)
			if (flags & chs::flags::build::REORDER)
			{
				const auto proj = [&](const auto & p) {
					const auto [i, j] = slice_.coord2indices(p);
					return slice_.indices2global(i, j);
				};
				const auto cmp = [&](const auto & a, const auto & b) { return proj(a) < proj(b); };
				if (flags & chs::flags::build::PARALLEL)
				{
					std::sort(std::execution::par_unseq, points.begin(), points.end(), cmp);
				}
				else { std::sort(points.begin(), points.end(), cmp); }
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

			const auto [min_i, min_j] = slice_.coord2indices(kernel.box().min());
			const auto [max_i, max_j] = slice_.coord2indices(kernel.box().max());

			for (const auto [i, j] :
			     ranges::views::cartesian_product(ranges::views::closed_indices(min_i, max_i),
			                                      ranges::views::closed_indices(min_j, max_j)))
			{
				const auto & cell_opt = slice_.at(i, j);
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
			const auto [p_i, p_j] = slice_.coord2indices(p);
			double search_radius  = slice_.idx2box(p_i, p_j).distance_to_wall(p, /* inside = */ true);

			auto candidates_within_sq_radius = [&] {
				const auto sq_radius = search_radius * search_radius;
				const auto it        = std::upper_bound(candidates.begin(), candidates.end(), sq_radius,
				                                        [](auto d, auto & pt) { return d < pt.first; });
				return it - candidates.begin();
			};

			// Taboo list (to avoid visiting the same cell twice)
			indices_array taboo_mins;
			indices_array taboo_maxs;

			auto is_tabooed = [&](const indices_array & indices) {
				return chs::within_closed_bounds<Dim>(indices, taboo_mins, taboo_maxs);
			};

			// Do an increasing search
			const double default_radius_increment = ranges::max(slice_.resolutions());

			// Explore first the neighbors of the cell containing p
			{
				taboo_mins = { p_i, p_j };
				taboo_maxs = { p_i, p_j };

				const auto & cell_opt = slice_.at(p_i, p_j);

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
				else { search_radius += default_radius_increment; }

				const auto [min_i, min_j] = slice_.coord2indices(p - search_radius);
				const auto [max_i, max_j] = slice_.coord2indices(p + search_radius);

				const indices_array min = { min_i, min_j };
				const indices_array max = { max_i, max_j };

				// If min == taboo_mins and max == taboo_maxs, we have already visited all the cells
				if (chs::all_equal<Dim>(min, taboo_mins) and chs::all_equal<Dim>(max, taboo_maxs))
				{
					continue;
				}

				for (const auto [i, j] :
				     ranges::views::cartesian_product(ranges::views::closed_indices(min_i, max_i),
				                                      ranges::views::closed_indices(min_j, max_j)))
				{
					if (is_tabooed({ i, j })) { continue; }

					const auto & cell_opt = slice_.at(i, j);

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
	};
} // namespace chs