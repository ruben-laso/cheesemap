#pragma once

#include <array>
#include <execution>
#include <queue>
#include <set>
#include <vector>

#include <range/v3/all.hpp>

#include "cheesemap/utils/Box.hpp"
#include "cheesemap/utils/Cell.hpp"

#include "cheesemap/concepts/concepts.hpp"

#include "cheesemap/kernels/kernels.hpp"

#include "cheesemap/maps/SmartSlice.hpp"

namespace chs
{
	template<typename Point_type>
	class Mixed3D
	{
		protected:
		static constexpr std::size_t Dim = 3;

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

		std::vector<chs::slice::Smart<Point_type>> slices_;

		[[nodiscard]] inline auto idx2box(const std::size_t i, const std::size_t j, const std::size_t k) const
		{
			Point center{
				box_.min()[0] + (static_cast<resolution_type>(i) + 0.5) * resolutions_[0], // x
				box_.min()[1] + (static_cast<resolution_type>(j) + 0.5) * resolutions_[1], // y
				box_.min()[2] + (static_cast<resolution_type>(k) + 0.5) * resolutions_[2]  // z
			};
			Point radii{ resolutions_[0] / 2, resolutions_[1] / 2, resolutions_[2] / 2 };

			return Box{ center, radii };
		}

		[[nodiscard]] inline auto coord2indices(const Point & p) const
		{
			const auto rel  = p - box_.min();
			const auto diff = box_.max() - box_.min();
			return std::make_tuple(
			        static_cast<std::size_t>(std::clamp(rel[0], 0.0, diff[0]) / resolutions_[0]),
			        static_cast<std::size_t>(std::clamp(rel[1], 0.0, diff[1]) / resolutions_[1]),
			        static_cast<std::size_t>(std::clamp(rel[2], 0.0, diff[2]) / resolutions_[2]));
		}

		[[nodiscard]] inline auto indices2global(const std::size_t i, const std::size_t j,
		                                         const std::size_t k) const
		{
			// return i * sizes_[1] * sizes_[2] + j * sizes_[2] + k;
			return k * sizes_[0] * sizes_[1] + j * sizes_[0] + i;
		}

		public:
		Mixed3D() = default;

		template<typename Points_rng>
		Mixed3D(Points_rng & points, const resolution_type res, const bool reorder = false) :
		        Mixed3D(points, dimensions_vector(Dim, res), reorder)
		{}

		template<typename Points_rng>
		Mixed3D(Points_rng & points, dimensions_vector res, const bool reorder = false) :
		        resolutions_(std::move(res)), box_(Box::mbb(points))
		{
			// Number of cells in each dimension
			const Point diff = box_.max() - box_.min();
			for (const auto i : ranges::views::indices(Dim))
			{
				sizes_[i] = static_cast<std::size_t>(std::floor(diff[i] / resolutions_[i])) + 1;
			}

			// Generate the slices
			for (const auto k : ranges::views::indices(sizes_[2]))
			{
				const auto box_min = idx2box(0, 0, k);
				const auto box_max = idx2box(sizes_[0] - 1, sizes_[1] - 1, k + 1);

				Box slice_box{ std::make_pair(box_min.min(), box_max.max()) };
				slices_.emplace_back(slice_box, resolutions_);
			}

			// Sort points by global idx (should improve locality when querying)
			if (reorder)
			{
				auto proj = [&](const auto & p) {
					const auto [i, j, k] = coord2indices(p);
					return indices2global(i, j, k);
				};
				std::sort(std::execution::par_unseq, points.begin(), points.end(),
				          [&](const auto & a, const auto & b) { return proj(a) < proj(b); });
			}

			// Add the points to the slices
			for (auto & point : points)
			{
				const auto [i, j, k] = coord2indices(point);
				slices_[k].add_point(point);
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

				for (const auto [i, j] :
				     ranges::views::cartesian_product(ranges::views::closed_indices(min_i, max_i),
				                                      ranges::views::closed_indices(min_j, max_j)))
				{
					const auto & cell_opt = slice.at(i, j);
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
			const auto [p_i, p_j, p_k] = coord2indices(p);
			double search_radius       = idx2box(p_i, p_j, p_k).closest_distance(p);

			auto candidates_within_radius = [&] {
				const auto it = std::upper_bound(candidates.begin(), candidates.end(), search_radius,
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
			const double default_radius_increment = ranges::max(resolutions_);

			// Explore first the neighbors of the cell containing p
			{
				taboo_mins = { p_i, p_j, p_k };
				taboo_maxs = { p_i, p_j, p_k };

				const auto & cell_opt = slices_[p_k].at(p_i, p_j);

				if (cell_opt.has_value())
				{
					ranges::for_each(cell_opt->get(), [&](const auto & point) {
						candidates.insert({ chs::distance(p, *point), point });
					});
				}
			}

			while (
			        // not enough candidates or last candidate is outside the search radius
			        (std::cmp_less(candidates.size(), k_neigh) or
			         candidates.back().first > search_radius) and
			        // we have not visited all the cells
			        not chs::all_visited<Dim>(taboo_mins, taboo_maxs, sizes_))
			{
				// Estimate the new required search radius -> k * density -> Saves time ~86% of the queries
				if (not candidates.empty() and search_radius > 0)
				{
					const auto density_based_radius = chs::radius_for_density(
					        candidates_within_radius(), search_radius, k_neigh);
					search_radius = std::min(density_based_radius,
					                         search_radius + default_radius_increment);
				}
				else { search_radius += default_radius_increment; }

				chs::kernels::Sphere<Dim> search(p, search_radius);

				const auto [min_i, min_j, min_k] = coord2indices(search.box().min());
				const auto [max_i, max_j, max_k] = coord2indices(search.box().max());

				for (const auto k : ranges::views::closed_indices(min_k, max_k))
				{
					const auto & slice = slices_[k];
					for (const auto [i, j] : ranges::views::cartesian_product(
					             ranges::views::closed_indices(min_i, max_i),
					             ranges::views::closed_indices(min_j, max_j)))
					{
						if (is_tabooed({ i, j, k })) { continue; }

						// Get the cell
						const auto & cell_opt = slice.at(i, j);

						if (not cell_opt.has_value()) { continue; }

						// If the cell is completely inside the search sphere, directly to candidates
						ranges::for_each(cell_opt->get(), [&](const auto & point) {
							candidates.insert({ chs::distance(p, *point), point });
						});
					}
				}

				taboo_mins = { min_i, min_j, min_k };
				taboo_maxs = { max_i, max_j, max_k };
			}

			return candidates;
		}
	};
} // namespace chs