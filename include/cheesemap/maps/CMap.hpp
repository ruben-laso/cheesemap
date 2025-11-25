#pragma once

#include <array>
#include <tuple>
#include <utility>
#include <vector>

#include <range/v3/all.hpp>

#include "cheesemap/utils/execution.hpp"

#include "cheesemap/concepts/concepts.hpp"

#include "cheesemap/utils/Box.hpp"
#include "cheesemap/utils/Cell.hpp"
#include "cheesemap/utils/sorted_vector.hpp"

#include "cheesemap/utils/arithmetic.hpp"
#include "cheesemap/utils/Cartesian.hpp"
#include "cheesemap/utils/flags.hpp"
#include "cheesemap/utils/type_traits.hpp"

#include "cheesemap/utils/inline.hpp"

namespace chs
{
	template<typename Point_type, std::size_t Dim, typename Derived_type>
	class CMap
	{
		protected:
		using resolution_type = double;
		using dimensions_type = chs::type_traits::tuple<resolution_type, Dim>;
		using gbl_idx_type    = std::size_t;
		using indices_type    = chs::type_traits::tuple<gbl_idx_type, Dim>;
		using cell_type       = Cell<Point_type>;

		static constexpr dimensions_type DEFAULT_RESOLUTIONS = chs::n_tuple<Dim>(resolution_type{ 1 });

		static constexpr double OVERALLOCATION_FACTOR = 1.2;

		// Dimension of each cell
		dimensions_type resolutions_ = DEFAULT_RESOLUTIONS;

		// Bounding box of the map
		Box box_;

		// Number of cells of the map on each dimension
		indices_type sizes_;

		// Cell (hyper-volume)
		double cell_volume_ = 1.0;

		// Density (pts/non-empty cell)
		double density_ = 1.0;

		template<gbl_idx_type... Is>
		[[nodiscard]] CHSINLINE auto indices2global(const auto & indices, std::index_sequence<Is...>) const
		{
			gbl_idx_type idx = 0;

			((idx = idx * std::get<Is>(sizes_) + std::get<Is>(indices)), ...);

			return idx;
		}

		[[nodiscard]] CHSINLINE auto indices2global(const auto & indices) const
		{
			return indices2global(indices, std::make_index_sequence<Dim>{});
		}

		template<gbl_idx_type... Is>
		[[nodiscard]] CHSINLINE auto idx2box(const auto & idx, std::index_sequence<Is...>) const
		{
			Point min = box_.min();
			Point max = box_.max();

			(((min[Is] = box_.min()[Is] +
			             static_cast<resolution_type>(std::get<Is>(idx)) * std::get<Is>(resolutions_)),
			  (max[Is] = min[Is] + std::get<Is>(resolutions_))),
			 ...);

			return Box(std::make_pair(min, max));
		}

		[[nodiscard]] CHSINLINE auto idx2box(const auto & idx) const
		{
			return idx2box(idx, std::make_index_sequence<Dim>{});
		}

		template<gbl_idx_type... Is>
		[[nodiscard]] CHSINLINE auto coord2indices(const Point & p, std::index_sequence<Is...>) const
		{
			indices_type idx;

			resolution_type diff;

			(((diff = p[Is] - box_.min()[Is]),
			  (diff = std::clamp(diff, 0.0, box_.max()[Is] - box_.min()[Is])),
			  (std::get<Is>(idx) = static_cast<gbl_idx_type>(diff / std::get<Is>(resolutions_)))),
			 ...);

			return idx;
		}

		[[nodiscard]] CHSINLINE auto coord2indices(const Point & p) const
		{
			return coord2indices(p, std::make_index_sequence<Dim>{});
		}

		[[nodiscard]] CHSINLINE auto estimate_pts_query(const Point_type &   kernel_min,
		                                                const Point_type &   kernel_max,
		                                                const indices_type & idx_min,
		                                                const indices_type & idx_max) const
		{
			// Estimate the required size of the output vector
			const auto query_bbox_vol = chs::volume_bbox<Dim>(kernel_min, kernel_max);
			const auto cells_bbox_vol =
			        static_cast<double>(chs::cartesian_product_size<Dim>(idx_min, idx_max)) * cell_volume_;
			const auto est_points = static_cast<std::size_t>((query_bbox_vol / cells_bbox_vol) * density_ *
			                                                 OVERALLOCATION_FACTOR);

			return est_points;
		}

		[[nodiscard]] auto cell_exists(const indices_type & indices) const -> bool
		{
			return static_cast<const Derived_type *>(this)->cell_exists_impl(indices);
		}

		[[nodiscard]] auto at(const indices_type & indices) -> cell_type &
		{
			return static_cast<Derived_type *>(this)->at_impl(indices);
		}

		[[nodiscard]] auto at(const indices_type & indices) const -> const cell_type &
		{
			return static_cast<const Derived_type *>(this)->at_impl(indices);
		}

		void add_point(const indices_type & idx, Point * point_ptr)
		{
			return static_cast<Derived_type *>(this)->add_point_impl(idx, point_ptr);
		}

		void allocate_cells() { return static_cast<Derived_type *>(this)->allocate_cells_impl(); }

		void shrink_to_fit() { return static_cast<Derived_type *>(this)->shrink_to_fit_impl(); }

		void init(std::vector<Point_type> & points, [[maybe_unused]] dimensions_type res,
		          chs::flags::build::flags_t flags = {})
		{
			cell_volume_ = [&]<gbl_idx_type... Is>(std::index_sequence<Is...>) {
				return (std::get<Is>(resolutions_) * ...);
			}(std::make_index_sequence<Dim>{});

			// Number of cells in each dimension
			[&]<gbl_idx_type... Is>(std::index_sequence<Is...>) {
				((std::get<Is>(sizes_) =
				          static_cast<gbl_idx_type>(std::floor((box_.max()[Is] - box_.min()[Is]) /
				                                               std::get<Is>(resolutions_))) +
				          1),
				 ...);
			}(std::make_index_sequence<Dim>{});

			// Create cells
			allocate_cells();

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

			// Assign points to cells
			ranges::for_each(points, [&](auto & point) { add_point(coord2indices(point), &point); });

			if (flags & chs::flags::build::SHRINK_TO_FIT) { shrink_to_fit(); }

			// Compute density
			const auto min             = chs::n_tuple<Dim>(0UL);
			const auto non_empty_cells = static_cast<std::size_t>(ranges::count_if(
			        chs::cartesian_open_bound<Dim>(min, sizes_), [&](const auto & indices) {
				        if (not cell_exists(indices)) { return false; }
				        const auto & cell = at(indices);
				        return not cell.empty();
			        }));
			density_ = static_cast<double>(points.size()) / static_cast<double>(non_empty_cells);
		}

		public:
		virtual ~CMap() = default;

		CMap() = delete;

		CMap(std::vector<Point_type> & points, const resolution_type res,
		     chs::flags::build::flags_t flags = {}) :
		        CMap(points, dimensions_type{ n_tuple<Dim>(res) }, flags)
		{}

		CMap(std::vector<Point_type> & points, dimensions_type res,
		     [[maybe_unused]] chs::flags::build::flags_t flags = {}) :
		        resolutions_(res), box_(Box::mbb(points))
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
			std::vector<Point *> points;

			const auto query_min = kernel.box().min();
			const auto query_max = kernel.box().max();

			const auto min = coord2indices(query_min);
			const auto max = coord2indices(query_max);

			// Estimate the required size of the output vector
			const auto est_points = estimate_pts_query(query_min, query_max, min, max);
			points.reserve(est_points);

			for (const auto & indices : chs::cartesian<Dim>(min, max))
			{
				if (not cell_exists(indices)) { continue; }

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

			// Taboo list (to avoid visiting the same cell twice)
			const auto center_idx = coord2indices(p);

			auto taboo_mins = chs::n_tuple<Dim>(std::numeric_limits<std::size_t>::max());
			auto taboo_maxs = chs::n_tuple<Dim>(std::numeric_limits<std::size_t>::max());

			auto is_taboo = [&](const auto & indices) {
				return chs::within_closed_bounds<Dim>(indices, taboo_mins, taboo_maxs);
			};

			// Do an increasing search
			const double default_radius_increment = chs::min<Dim>(resolutions_);

			// Search radius starts within the cell containing p
			auto search_radius = 0.0;
			auto next_radius   = idx2box(coord2indices(p)).distance_to_wall(p, /* inside = */ true);

			auto candidates_within_sq_radius = [&] {
				const auto sq_radius = search_radius * search_radius;
				const auto it        = std::upper_bound(candidates.begin(), candidates.end(), sq_radius,
				                                        [](auto d, auto & pt) { return d < pt.first; });
				return it - candidates.begin();
			};

			while (
			        // not enough candidates or last candidate is outside the search radius
			        (std::cmp_less(candidates.size(), k) or
			         candidates.back().first > (search_radius * search_radius)) and
			        // we have not visited all the cells
			        not chs::all_visited<Dim>(taboo_mins, taboo_maxs, sizes_))
			{
				search_radius = next_radius;

				const auto min = coord2indices(p - search_radius);
				const auto max = coord2indices(p + search_radius);

				// If min == taboo_mins and max == taboo_maxs, we have already visited all those cells
				// If not, we search those cells
				if (not(chs::all_equal<Dim>(min, taboo_mins) and chs::all_equal<Dim>(max, taboo_maxs)))
				{
					for (const auto & indices : chs::cartesian<Dim>(min, max))
					{
						if (is_taboo(indices)) { continue; }
						if (not cell_exists(indices)) { continue; }

						const auto & cell = at(indices);

						for (const auto & point : cell)
						{
							candidates.insert({ chs::sq_distance(p, *point), point });
						}
					}
				}

				// Update the taboo indices
				taboo_mins = min;
				taboo_maxs = max;

				// Estimate the new required search radius -> k * density -> Saves time ~86% of the queries
				if (not candidates.empty() and search_radius > 0)
				{
					const auto density_based_radius = chs::radius_for_density(
					        candidates_within_sq_radius(), search_radius, k);
					next_radius = std::min(density_based_radius,
					                       search_radius + default_radius_increment);
				}
				else { next_radius = search_radius + default_radius_increment; }
			}

			ranges::for_each(candidates,
			                 [](auto & candidate) { candidate.first = std::sqrt(candidate.first); });

			return candidates;
		}

		// [[nodiscard]] inline auto cells_stored() const { return std::vector<bool>(cells_.size(), true); }
		//
		// [[nodiscard]] inline auto points_per_cell() const
		// {
		// 	std::vector<std::size_t> num_points(cells_.size());
		// 	ranges::transform(cells_, num_points.begin(), [](const auto & cell) { return cell.size(); });
		// 	return num_points;
		// }
		//
		// [[nodiscard]] inline auto get_num_cells() const { return cells_.size(); }
		//
		// [[nodiscard]] inline auto get_num_empty_cells() const
		// {
		// 	return ranges::count_if(cells_, [](const auto & cell) { return cell.empty(); });
		// }
	};
} // namespace chs