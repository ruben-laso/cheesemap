#pragma once

#include <array>
#include <optional>
#include <vector>

#include "cheesemap/utils/Box.hpp"
#include "cheesemap/utils/Cell.hpp"
#include "cheesemap/utils/type_traits.hpp"

namespace chs::slice
{
	template<typename Point_type, template<typename, typename...> class HashMap = std::unordered_map>
	class Smart
	{
		protected:
		static constexpr auto SPARSE_TO_DENSE_THRESHOLD = 0.8;

		static constexpr auto Dim = 2;

		using resolution_type = double;
		using dimensions_type = chs::type_traits::tuple<resolution_type, Dim>;
		using indices_type    = chs::type_traits::tuple<std::size_t, Dim>;
		using cell_type       = Cell<Point_type>;

		static constexpr dimensions_type DEFAULT_RESOLUTIONS = chs::n_tuple<Dim>(resolution_type{ 1 });

		// Dimension of each cell
		dimensions_type resolutions_ = DEFAULT_RESOLUTIONS;

		// Bounding box of the map
		Box box_;

		// Number of cells of the map on each dimension
		indices_type sizes_;

		// Cells of the map (using sparse representation)
		HashMap<std::size_t, cell_type> cells_sparse_;

		// Cells of the map (using dense representation)
		std::vector<cell_type> cells_dense_;

		// Use sparse
		bool use_sparse_ = true;

		inline void add_point_dense(Point_type & point)
		{
			const auto [i, j] = coord2indices(point);
			const auto idx    = i * std::get<1>(sizes_) + j;
			cells_dense_[idx].emplace_back(&point);
		}

		inline void add_point_sparse(Point_type & point)
		{
			const auto idx = indices2global(coord2indices(point));

			auto cells_it = cells_sparse_.find(idx);

			if (cells_it == cells_sparse_.end())
			{
				cells_sparse_.emplace(idx, cell_type{});
				cells_it = cells_sparse_.find(idx);
			}

			cells_it->second.emplace_back(&point);

			if (density() > SPARSE_TO_DENSE_THRESHOLD) { sparse2dense(); }
		}

		inline void sparse2dense()
		{
			cells_dense_ = {};

			for (const auto indices :
			     ranges::views::cartesian_product(ranges::views::indices(std::get<0>(sizes_)),
			                                      ranges::views::indices(std::get<1>(sizes_))))
			{
				const auto idx = indices2global(indices);
				if (const auto cell_it = cells_sparse_.find(idx); cell_it != cells_sparse_.end())
				{
					cells_dense_.emplace_back(std::move(cell_it->second));
				}
				else { cells_dense_.emplace_back(cell_type{}); }
			}

			use_sparse_   = false;
			cells_sparse_ = {};
		}

		inline void dense2sparse()
		{
			cells_sparse_ = {};

			for (const auto indices :
			     ranges::views::cartesian_product(ranges::views::indices(std::get<0>(sizes_)),
			                                      ranges::views::indices(std::get<1>(sizes_))))
			{
				const auto idx = indices2global(indices);
				if (not cells_dense_[idx].empty())
				{
					cells_sparse_.emplace(idx, std::move(cells_dense_[idx]));
				}
			}

			use_sparse_ = true;

			cells_dense_ = {};
		}

		public:
		Smart() = delete;

		Smart(const Box & box, const resolution_type res) : Smart(box, dimensions_type(chs::n_tuple<Dim>(res)))
		{}

		Smart(const Box & box, const dimensions_type & res) :
		        resolutions_(res),
		        box_(box),
		        sizes_({ static_cast<std::size_t>(std::ceil((box.max()[0] - box.min()[0]) / std::get<0>(res))),
		                 static_cast<std::size_t>(std::ceil((box.max()[1] - box.min()[1]) / std::get<1>(res))) })
		{}

		template<typename Points_rng>
		        requires ranges::input_range<Points_rng>
		Smart(Points_rng & points, const resolution_type res) : Smart(Box::mbb(points), res)
		{
			for (auto & point : points)
			{
				add_point(point);
			}
		}

		[[nodiscard]] inline auto sparse() const { return use_sparse_; }

		[[nodiscard]] inline auto resolutions() const -> const auto & { return resolutions_; }

		[[nodiscard]] inline auto resolutions() { return resolutions_; }

		[[nodiscard]] inline auto sizes() const -> const auto & { return sizes_; }

		[[nodiscard]] inline auto size() const { return std::get<0>(sizes_) * std::get<1>(sizes_); }

		[[nodiscard]] inline auto density() const
		{
			if (use_sparse_)
			{
				return static_cast<double>(cells_sparse_.size()) / static_cast<double>(size());
			}
			return 1.0;
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

		[[nodiscard]] inline auto indices2global(const indices_type & indices) const
		{
			return std::get<0>(indices) * std::get<1>(sizes_) + std::get<1>(indices);
		}

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

		inline void add_point(Point_type & point)
		{
			if (use_sparse_) { add_point_sparse(point); }
			else { add_point_dense(point); }
		}

		inline void shrink_to_fit()
		{
			if (use_sparse_)
			{
				ranges::for_each(cells_sparse_, [](auto & cell) { cell.second.shrink_to_fit(); });
			}
			else
			{
				ranges::for_each(cells_dense_, [](auto & cell) { cell.shrink_to_fit(); });
			}
		}

		[[nodiscard]] CHSINLINE auto find(const indices_type & indices) -> cell_type *
		{
			const auto & gbl_idx = indices2global(indices);

			if (use_sparse_)
			{
				const auto it = cells_sparse_.find(gbl_idx);
				return it == cells_sparse_.end() ? nullptr : &(it->second);
			}
			return &cells_dense_[gbl_idx];
		}

		[[nodiscard]] CHSINLINE auto find(const indices_type & indices) const -> const cell_type *
		{
			const auto gbl_idx = indices2global(indices);
			if (use_sparse_)
			{
				const auto it = cells_sparse_.find(gbl_idx);
				return it == cells_sparse_.end() ? nullptr : &(it->second);
			}
			return &cells_dense_[gbl_idx];
		}

		[[nodiscard]] CHSINLINE auto at(const indices_type & indices) -> cell_type * { return find(indices); }

		[[nodiscard]] CHSINLINE auto at(const indices_type & indices) const -> const cell_type *
		{
			return find(indices);
		}
	};
} // namespace chs::slice