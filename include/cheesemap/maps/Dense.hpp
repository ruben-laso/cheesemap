#pragma once

#include "cheesemap/maps/IMap.hpp"

namespace chs
{
	template<typename Point_type, std::size_t Dim = 3>
	class Dense : public IMap<Point_type, Dim>
	{
		using resolution_type = double;
		using dimensions_type = chs::type_traits::tuple<resolution_type, Dim>;
		using gbl_idx_type    = std::size_t;
		using indices_type    = chs::type_traits::tuple<gbl_idx_type, Dim>;
		using cell_type       = Cell<Point_type>;

		// Cells of the map
		std::vector<cell_type> cells_;

		[[nodiscard]] inline auto cell_exists([[maybe_unused]] const indices_type & indices) const
		        -> bool override
		{
			return true;
		}

		[[nodiscard]] inline auto at(const indices_type & indices) -> cell_type & override
		{
			return cells_[this->indices2global(indices)];
		}

		[[nodiscard]] inline auto at(const indices_type & indices) const -> const cell_type & override
		{
			return cells_[this->indices2global(indices)];
		}

		void add_point(const indices_type & idx, Point * point_ptr) override
		{
			cells_[this->indices2global(idx)].emplace_back(point_ptr);
		}

		inline void allocate_cells() override
		{
			const auto num_cells = [&]<std::size_t... Is>(std::index_sequence<Is...>) {
				std::size_t acc = 1;
				((acc *= std::get<Is>(Dense::sizes_)), ...);
				return acc;
			}(std::make_index_sequence<Dim>{});

			cells_.resize(num_cells);
		}

		void shrink_to_fit() override
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
		        IMap<Point_type, Dim>(points, res, flags)
		{
			this->init(points, res, flags);
		}

		[[nodiscard]] inline auto cells_stored() const { return std::vector<bool>(cells_.size(), true); }

		[[nodiscard]] inline auto points_per_cell() const
		{
			std::vector<std::size_t> num_points(cells_.size());
			ranges::transform(cells_, num_points.begin(), [](const auto & cell) { return cell.size(); });
			return num_points;
		}

		[[nodiscard]] inline auto get_num_cells() const { return cells_.size(); }

		[[nodiscard]] inline auto get_num_empty_cells() const
		{
			return ranges::count_if(cells_, [](const auto & cell) { return cell.empty(); });
		}
	};
} // namespace chs