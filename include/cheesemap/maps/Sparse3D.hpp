#pragma once

#include <array>
#include <queue>
#include <set>
#include <vector>

#include "cheesemap/maps/Sparse.hpp"

namespace chs
{
	template<typename Point_type>
	class Sparse3D : public Sparse<Point_type, 3>
	{
		protected:
		static constexpr std::size_t Dim = 3;

		using super_type = chs::Sparse<Point_type, Dim>;

		public:
		Sparse3D() = delete;

		template<typename... Args_types>
		Sparse3D(Args_types &&... args) : super_type(std::forward<Args_types>(args)...)
		{}

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

			const auto min = this->coord2indices(kernel.box().min());
			const auto max = this->coord2indices(kernel.box().max());

			for (const auto [i, j, k] :
			     ranges::views::cartesian_product(ranges::views::closed_indices(min[0], max[0]),
			                                      ranges::views::closed_indices(min[1], max[1]),
			                                      ranges::views::closed_indices(min[2], max[2])))
			{
				const auto global_idx = i * this->sizes_[1] * this->sizes_[2] + j * this->sizes_[2] + k;

				const auto cell_it = this->cells_.find(global_idx);

				if (cell_it == this->cells_.end()) { continue; }

				const auto & cell = cell_it->second;

				ranges::for_each(cell.points(), [&](auto * point_ptr) {
					if (kernel.is_inside(*point_ptr) && filter(*point_ptr))
					{
						points.emplace_back(point_ptr);
					}
				});
			}

			return points;
		}

		[[nodiscard]] inline auto knn(const std::integral auto k_neigh, const Point_type & p) const
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
			double search_radius = ranges::max(this->resolutions_);

			const auto num_cells =
			        ranges::accumulate(this->sizes_, std::size_t{ 1 }, std::multiplies<std::size_t>{});

			while (std::cmp_less(candidates.size(), k_neigh) and std::cmp_less(taboo.size(), num_cells))
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
				const auto min = this->coord2indices(search.box().min());
				const auto max = this->coord2indices(search.box().max());

				for (const auto [i, j, k] :
				     ranges::views::cartesian_product(ranges::views::closed_indices(min[0], max[0]),
				                                      ranges::views::closed_indices(min[1], max[1]),
				                                      ranges::views::closed_indices(min[2], max[2])))
				{
					const auto global_idx =
					        i * this->sizes_[1] * this->sizes_[2] + j * this->sizes_[2] + k;

					// Skip already visited cells
					if (taboo.contains(global_idx)) { continue; }

					// Mark the cell as visited
					taboo.insert(global_idx);

					// Get the cell
					const auto cell_it = this->cells_.find(global_idx);
					if (cell_it == this->cells_.end()) { continue; }

					const auto & cell = cell_it->second;

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

			return candidates | ranges::views::take(k_neigh) | ranges::to_vector;
		}
	};
} // namespace chs