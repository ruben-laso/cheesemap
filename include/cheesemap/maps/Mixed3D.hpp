#pragma once

#include <unordered_map>
#include <vector>

#include "cheesemap/maps/CMap.hpp"

#include "cheesemap/utils/inline.hpp"

namespace chs
{
	template<typename Point_type, template<typename, typename...> class HashMap = std::unordered_map>
	class Mixed3D : public CMap<Point_type, 3, Mixed3D<Point_type, HashMap>>
	{
		protected:
		static constexpr std::size_t Dim = 3;

		using Super = CMap<Point_type, 3, Mixed3D>;
		friend class CMap<Point_type, 3, Mixed3D>;

		using resolution_type = typename Super::resolution_type;
		using dimensions_type = typename Super::dimensions_type;
		using indices_type    = typename Super::indices_type;
		using cell_type       = typename Super::cell_type;

		std::vector<chs::slice::Smart<Point_type, HashMap>> slices_;

		[[nodiscard]] CHSINLINE auto find(const indices_type & indices) const -> const cell_type *
		{
			const auto & [i, j, k] = indices;
			return slices_[k].find({ i, j });
		}

		[[nodiscard]] CHSINLINE auto find(const indices_type & indices) -> cell_type *
		{
			const auto & [i, j, k] = indices;
			return slices_[k].find({ i, j });
		}

		[[nodiscard]] CHSINLINE auto cell_exists_impl(const indices_type & indices) const -> bool
		{
			return find(indices) != nullptr;
		}

		[[nodiscard]] CHSINLINE auto at_impl(const indices_type & indices) -> cell_type &
		{
			return *find(indices);
		}

		[[nodiscard]] CHSINLINE auto at_impl(const indices_type & indices) const -> const cell_type &
		{
			return *find(indices);
		}

		CHSINLINE void add_point_impl(const indices_type & idx, Point * point_ptr)
		{
			const auto & [i, j, k] = idx;
			slices_[k].add_point(*point_ptr);
		}

		void allocate_cells_impl()
		{
			// Generate the slices
			for (const auto k : ranges::views::indices(std::get<2>(this->sizes_)))
			{
				const auto box_min = this->idx2box(indices_type{ 0, 0, k });
				const auto box_max = this->idx2box(indices_type(std::get<0>(this->sizes_) - 1,
				                                                std::get<1>(this->sizes_) - 1, k + 1));

				Box slice_box{ std::make_pair(box_min.min(), box_max.max()) };
				slices_.emplace_back(slice_box, std::make_tuple(std::get<0>(this->resolutions_),
				                                                std::get<1>(this->resolutions_)));
			}
		}

		void shrink_to_fit_impl()
		{
			for (auto & slice : slices_)
			{
				slice.shrink_to_fit();
			}
		}

		public:
		Mixed3D() = default;

		template<typename Points_rng>
		Mixed3D(Points_rng & points, const resolution_type res, const chs::flags::build::flags_t flags = {}) :
		        Mixed3D(points, dimensions_type(n_tuple<Dim>(res)), flags)
		{}

		template<typename Points_rng>
		Mixed3D(Points_rng & points, dimensions_type res, const chs::flags::build::flags_t flags = {}) :
		        CMap<Point_type, 3, Mixed3D>(points, res, flags)
		{
			this->init(points, res, flags);
		}
	};
} // namespace chs