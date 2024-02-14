#pragma once

#include <array>
#include <queue>
#include <set>
#include <vector>

#include <armadillo>

#include "cheesemap/concepts/concepts.hpp"

#include "cheesemap/kernels/Sphere.hpp"
#include "cheesemap/utils/Cell.hpp"

namespace chs
{
	template<typename Point_type>
	class Sparse3D
	{
		protected:
		static constexpr std::size_t Dim = 3;

		using resolution_type   = double;
		using dimensions_array  = std::array<resolution_type, Dim>;
		using dimensions_vector = std::vector<resolution_type>;
		using indices_array     = std::array<std::size_t, Dim>;
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

		// Sparse matrix storing the indices of the cells
		arma::Cube<std::size_t> indices_;

		// Cells of the map
		std::vector<cell_type> cells_;

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
			const auto size = box_.max() - box_.min();
			const auto i    = static_cast<std::size_t>(std::clamp(rel[0], 0.0, size[0]) / resolutions_[0]);
			const auto j    = static_cast<std::size_t>(std::clamp(rel[1], 0.0, size[1]) / resolutions_[1]);
			const auto k    = static_cast<std::size_t>(std::clamp(rel[2], 0.0, size[2]) / resolutions_[2]);
			return std::make_tuple(i, j, k);
		}

		[[nodiscard]] inline auto submat(const Point_type & min, const Point_type & max) const
		{
			const auto [min_i, min_j, min_k] = coord2indices(min);
			const auto [max_i, max_j, max_k] = coord2indices(max);

			return std::as_const(indices_).subcube(min_i, min_j, min_k, max_i, max_j, max_k);
		}

		[[nodiscard]] auto global_idx_cells_to_search(auto && kernel) const
		{
			std::vector<std::size_t> cells_to_search;

			auto indices = submat(kernel.box().min(), kernel.box().max());

			return indices | ranges::views::filter([](const auto & global_idx) {
				       return std::cmp_greater(global_idx, 0);
			       }) |
			       ranges::views::transform([](const auto & global_idx) { return global_idx - 1; }) |
			       ranges::to_vector;
		}

		public:
		Sparse3D() = delete;

		template<typename Points_rng>
		Sparse3D(Points_rng & points, const resolution_type res) : Sparse3D(points, dimensions_vector(Dim, res))
		{}

		template<typename Points_rng>
		Sparse3D(Points_rng & points, dimensions_vector res) :
		        resolutions_(std::move(res)), box_(Box::mbb(points))
		{
			for (const auto i : ranges::views::indices(Dim))
			{
				sizes_[i] = static_cast<std::size_t>(
				        std::floor((box_.max()[i] - box_.min()[i]) / resolutions_[i]) + 1);
			}

			indices_.resize(sizes_[0], sizes_[1], sizes_[2]);

			for (auto & point : points)
			{
				const auto [i, j, k] = coord2indices(point);

				// If the cell is empty
				if (std::as_const(indices_).at(i, j, k) == 0)
				{
					// Create a new cell
					cells_.emplace_back(Cell<Point_type>{ idx2box(i, j, k) });
					// Store the index of the cell in the sparse matrix
					indices_.at(i, j, k) = cells_.size();
				}

				// Insert in (i, j) the index of the cell - 1 (0 is reserved for empty cells)
				const auto idx = std::as_const(indices_).at(i, j, k) - 1;
				cells_[idx].add_point(&point);
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

			const auto & global_indices = submat(kernel.box().min(), kernel.box().max());

			for (const auto global_idx : global_indices)
			{
				// If 0 -> empty cell
				if (std::cmp_equal(global_idx, 0)) { continue; }

				const auto & cell = cells_[global_idx - 1];

				ranges::for_each(cell.points(), [&](auto * point_ptr) {
					if (kernel.is_inside(*point_ptr) && filter(*point_ptr))
					{
						points.emplace_back(point_ptr);
					}
				});
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

			auto not_visited = [&](const auto & global_idx) { return not taboo.contains(global_idx); };

			// Do an increasing search
			double search_radius = ranges::max(resolutions_);

			while (std::cmp_less(candidates.size(), k) and std::cmp_less(taboo.size(), cells_.size()))
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
				const auto cells_to_search = global_idx_cells_to_search(search);

				// Visit the (non-visited already) cells
				for (const auto cell_idx : cells_to_search | ranges::views::filter(not_visited))
				{
					// Mark the cell as visited
					taboo.insert(cell_idx);

					// Get the cell
					auto & cell = cells_[cell_idx];

					// If the cell is completely inside the search sphere, directly to candidates
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

			return candidates | ranges::views::take(k) | ranges::to_vector;
		}
	};
} // namespace chs