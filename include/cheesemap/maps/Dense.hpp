#pragma once

#include <vector>

#include "cheesemap/maps/CMap.hpp"

#include "cheesemap/utils/inline.hpp"

namespace chs
{
	template<typename Point_type, std::size_t Dim = 3>
	class Dense : public CMap<Point_type, Dim, Dense<Point_type, Dim>>
	{
		using Super = CMap<Point_type, Dim, Dense>;
		// grant CMap access
		friend class CMap<Point_type, Dim, Dense>;

		using resolution_type = typename Super::resolution_type;
		using dimensions_type = typename Super::dimensions_type;
		using gbl_idx_type    = typename Super::gbl_idx_type;
		using indices_type    = typename Super::indices_type;
		using cell_type       = typename Super::cell_type;

		// Cells of the map
		std::vector<cell_type> cells_;

		[[nodiscard]] CHSINLINE auto cell_exists_impl([[maybe_unused]] const indices_type & indices) const
		        -> bool
		{
			return true;
		}

		[[nodiscard]] CHSINLINE auto at_impl(const indices_type & indices) -> cell_type &
		{
			return cells_[this->indices2global(indices)];
		}

		[[nodiscard]] CHSINLINE auto at_impl(const indices_type & indices) const -> const cell_type &
		{
			return cells_[this->indices2global(indices)];
		}

		CHSINLINE void add_point_impl(const indices_type & idx, Point * point_ptr)
		{
			cells_[this->indices2global(idx)].emplace_back(point_ptr);
		}

		void allocate_cells_impl()
		{
			const auto num_cells = [&]<std::size_t... Is>(std::index_sequence<Is...>) {
				std::size_t acc = 1;
				((acc *= std::get<Is>(Super::sizes_)), ...);
				return acc;
			}(std::make_index_sequence<Dim>{});

			cells_.resize(num_cells);
		}

		void shrink_to_fit_impl()
		{
			for (auto & cell : cells_)
			{
				cell.shrink_to_fit();
			}
		}

		public:
		Dense() = delete;

		Dense(std::vector<Point_type> & points, const resolution_type res,
		      chs::flags::build::flags_t flags = {}) :
		        Dense(points, dimensions_type{ n_tuple<Dim>(res) }, flags)
		{}

		Dense(std::vector<Point_type> & points, dimensions_type res, chs::flags::build::flags_t flags = {}) :
		        Super(points, res, flags)
		{
			this->init(points, res, flags);
		}

		[[nodiscard]] CHSINLINE auto cells_stored() const { return std::vector<bool>(cells_.size(), true); }

		[[nodiscard]] CHSINLINE auto points_per_cell() const
		{
			std::vector<std::size_t> num_points(cells_.size());
			ranges::transform(cells_, num_points.begin(), [](const auto & cell) { return cell.size(); });
			return num_points;
		}

		[[nodiscard]] CHSINLINE auto get_num_cells() const { return cells_.size(); }

		[[nodiscard]] CHSINLINE auto get_num_empty_cells() const
		{
			return ranges::count_if(cells_, [](const auto & cell) { return cell.empty(); });
		}
	};
} // namespace chs