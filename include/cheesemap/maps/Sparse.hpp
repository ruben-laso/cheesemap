#pragma once

#include <unordered_map>

#include "cheesemap/maps/CMap.hpp"

namespace chs
{
	template<typename Point_type, std::size_t Dim = 3,
	         template<typename, typename...> class HashMap = std::unordered_map>
	class Sparse : public CMap<Point_type, Dim, Sparse<Point_type, Dim, HashMap>>
	{
		using Super = CMap<Point_type, Dim, Sparse<Point_type, Dim, HashMap>>;
		friend class CMap<Point_type, Dim, Sparse<Point_type, Dim, HashMap>>;

		protected:
		using resolution_type = typename Super::resolution_type;
		using dimensions_type = typename Super::dimensions_type;
		using indices_type    = typename Super::indices_type;
		using cell_type       = typename Super::cell_type;

		// Cells of the map
		HashMap<std::size_t, cell_type> cells_;

		[[nodiscard]] CHSINLINE auto find(const indices_type & indices)
		{
			return cells_.find(this->indices2global(indices));
		}

		[[nodiscard]] CHSINLINE auto find(const indices_type & indices) const
		{
			return cells_.find(this->indices2global(indices));
		}

		[[nodiscard]] CHSINLINE auto cell_exists_impl(const indices_type & indices) const -> bool
		{
			return find(indices) != cells_.end();
		}

		[[nodiscard]] CHSINLINE auto at_impl(const indices_type & indices) -> cell_type &
		{
			return find(indices)->second;
		}

		[[nodiscard]] CHSINLINE auto at_impl(const indices_type & indices) const -> const cell_type &
		{
			return find(indices)->second;
		}

		CHSINLINE void add_point_impl(const indices_type & idx, Point * point_ptr)
		{
			cells_[this->indices2global(idx)].emplace_back(point_ptr);
		}

		void allocate_cells_impl()
		{
			// Nothing to do
		}

		void shrink_to_fit_impl()
		{
			for (auto & [idx, cell] : cells_)
			{
				cell.shrink_to_fit();
			}
		}

		public:
		Sparse() = delete;

		template<typename Points_rng>
		Sparse(Points_rng & points, const resolution_type res, const chs::flags::build::flags_t flags = {}) :
		        Sparse(points, dimensions_type(n_tuple<Dim>(res)), flags)
		{}

		template<typename Points_rng>
		Sparse(Points_rng & points, dimensions_type res, const chs::flags::build::flags_t flags = {}) :
		        Super(points, res, flags)
		{
			this->init(points, res, flags);
		}

		[[nodiscard]] CHSINLINE auto cells_stored() const
		{
			std::vector<bool> cells_stored(chs::product<Dim>(this->sizes_), false);
			for (const auto & [idx, cell] : cells_)
			{
				cells_stored[idx] = true;
			}
			return cells_stored;
		}

		[[nodiscard]] CHSINLINE auto points_per_cell() const
		{
			std::vector<std::size_t> num_points(chs::product<Dim>(this->sizes_));
			for (const auto & [idx, cell] : cells_)
			{
				num_points[idx] = cell.size();
			}
			return num_points;
		}

		[[nodiscard]] CHSINLINE auto get_num_cells() const { return cells_.size(); }

		[[nodiscard]] CHSINLINE auto get_num_empty_cells() const
		{
			return ranges::count_if(cells_, [](const auto & cell) { return cell.second.empty(); });
		}
	};
} // namespace chs