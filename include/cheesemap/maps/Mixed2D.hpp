#pragma once

#include <execution>
#include <set>
#include <vector>

#include <range/v3/all.hpp>

#include "cheesemap/concepts/concepts.hpp"

#include "cheesemap/kernels/kernels.hpp"

#include "cheesemap/maps/SmartSlice.hpp"

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
		Mixed2D(Points_rng & points, const resolution_type res, const bool reorder = false) :
		        Mixed2D(points, dimensions_vector(Dim, res), reorder)
		{}

		template<typename Points_rng>
		Mixed2D(Points_rng & points, const dimensions_vector & res, const bool reorder = false) :
		        slice_(Box::mbb(points), res)
		{
			// Sort points by global idx (should improve locality when querying)
			if (reorder)
			{
				auto proj = [&](const auto & p) {
					const auto [i, j] = slice_.coord2indices(p);
					return slice_.indices2global(i, j);
				};
#ifdef CHS_HAVE_EXECUTION
				std::sort(std::execution::par_unseq, points.begin(), points.end(),
				          [&](const auto & a, const auto & b) { return proj(a) < proj(b); });
#else
				std::sort(points.begin(), points.end(),
				          [&](const auto & a, const auto & b) { return proj(a) < proj(b); });
#endif
			}

			for (auto & point : points)
			{
				slice_.add_point(point);
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
			const auto [p_i, p_j] = slice_.coord2indices(p);
			double search_radius  = slice_.idx2box(p_i, p_j).closest_distance(p);

			auto candidates_within_radius = [&] {
				const auto it = std::upper_bound(candidates.begin(), candidates.end(), search_radius,
				                                 [](auto d, auto & pt) { return d < pt.first; });
				return it - candidates.begin();
			};

			// Taboo list (to avoid visiting the same cell twice)
			indices_array taboo_mins;
			indices_array taboo_maxs;

			auto is_tabooed = [&](const auto i, const auto j) {
				auto within_bounds = [](const auto idx, const auto min, const auto max) {
					return std::cmp_greater_equal(idx, min) and std::cmp_less_equal(idx, max);
				};
				return within_bounds(i, taboo_mins[0], taboo_maxs[0]) and
				       within_bounds(j, taboo_mins[1], taboo_maxs[1]);
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
					ranges::for_each(cell_opt->get(), [&](const auto & point_ptr) {
						const auto d = distance(p, *point_ptr);
						candidates.insert({ d, point_ptr });
					});
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
			        ranges::any_of(ranges::views::zip(taboo_mins, taboo_maxs, slice_.sizes()),
			                       [](const auto & t) {
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

				const auto [min_i, min_j] = slice_.coord2indices(search.box().min());
				const auto [max_i, max_j] = slice_.coord2indices(search.box().max());

				for (const auto [i, j] :
				     ranges::views::cartesian_product(ranges::views::closed_indices(min_i, max_i),
				                                      ranges::views::closed_indices(min_j, max_j)))
				{
					if (is_tabooed(i, j)) { continue; }

					// Get the cell
					const auto & cell_opt = slice_.at(i, j);

					if (not cell_opt.has_value()) { continue; }

					// If the cell is completely inside the search sphere, directly to candidates
					ranges::for_each(cell_opt->get(), [&](const auto & point_ptr) {
						const auto d = distance(p, *point_ptr);
						candidates.insert({ d, point_ptr });
					});
				}

				taboo_mins = { min_i, min_j };
				taboo_maxs = { max_i, max_j };
			}

			return candidates;
		}
	};
} // namespace chs