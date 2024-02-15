#pragma once

#include <queue>
#include <set>
#include <vector>

#include <range/v3/all.hpp>

#include "cheesemap/concepts/concepts.hpp"
#include "cheesemap/kernels/kernels.hpp"
#include "cheesemap/maps/SmartSlice.hpp"

namespace chs
{
	template<typename Point_type>
	class Mixed2D
	{
		protected:
		static constexpr std::size_t Dim = 2;

		using resolution_type   = double;
		using dimensions_vector = std::vector<resolution_type>;

		Slice::Smart<Point_type> slice_;

		public:
		Mixed2D() = default;

		template<typename Points_rng>
		Mixed2D(Points_rng & points, const resolution_type res) : Mixed2D(points, dimensions_vector(Dim, res))
		{}

		template<typename Points_rng>
		Mixed2D(Points_rng & points, const dimensions_vector & res) : slice_(Box::mbb(points), res)
		{
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
				for (auto * point_ptr : cell_opt->get().points())
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

			auto is_visited = [&](const auto i, const auto j) {
				const auto global_idx = slice_.to_global_idx_as_dense(i, j);
				return taboo.contains(global_idx);
			};

			// Do an increasing search
			double search_radius = ranges::max(slice_.resolutions());

			while (std::cmp_less(candidates.size(), k) and std::cmp_less(taboo.size(), slice_.size()))
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

				const auto [min_i, min_j] = slice_.coord2indices(search.box().min());
				const auto [max_i, max_j] = slice_.coord2indices(search.box().max());

				for (const auto [i, j] :
				     ranges::views::cartesian_product(ranges::views::closed_indices(min_i, max_i),
				                                      ranges::views::closed_indices(min_j, max_j)))
				{
					// Slip already visited cells
					if (is_visited(i, j)) { continue; }

					// Mark the cell as visited
					taboo.insert(slice_.to_global_idx_as_dense(i, j));

					// Get the cell
					const auto & cell_opt = slice_.at(i, j);

					if (not cell_opt.has_value()) { continue; }

					// If the cell is completely inside the search sphere, directly to candidates
					ranges::for_each(cell_opt->get().points(), [&](const auto & point_ptr) {
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