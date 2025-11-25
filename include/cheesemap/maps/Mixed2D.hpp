#pragma once

#include <unordered_map>

#include "cheesemap/maps/SmartSlice.hpp"

#include "cheesemap/maps/CMap.hpp"

#include "cheesemap/utils/inline.hpp"

namespace chs
{
	template<typename Point_type, template<typename, typename...> class HashMap = std::unordered_map>
	class Mixed2D : public CMap<Point_type, 2, Mixed2D<Point_type, HashMap>>
	{
		protected:
		static constexpr std::size_t Dim = 2;

		using Super = CMap<Point_type, Dim, Mixed2D>;
		friend class CMap<Point_type, Dim, Mixed2D>;

		using resolution_type = typename Super::resolution_type;
		using dimensions_type = typename Super::dimensions_type;
		using indices_type    = typename Super::indices_type;
		using cell_type       = typename Super::cell_type;

		chs::slice::Smart<Point_type, HashMap> slice_;

		[[nodiscard]] CHSINLINE auto find(const indices_type & indices) const { return slice_.find(indices); }

		[[nodiscard]] CHSINLINE auto cell_exists_impl(const indices_type & indices) const -> bool
		{
			return slice_.find(indices) != nullptr;
		}

		[[nodiscard]] CHSINLINE auto at_impl(const indices_type & indices) -> cell_type &
		{
			return *slice_.at(indices);
		}

		[[nodiscard]] CHSINLINE auto at_impl(const indices_type & indices) const -> const cell_type &
		{
			return *slice_.at(indices);
		}

		CHSINLINE void add_point_impl([[maybe_unused]] const indices_type & idx, Point * point_ptr)
		{
			slice_.add_point(*point_ptr);
		}

		void allocate_cells_impl()
		{
			// Nothing to do
		}

		void shrink_to_fit_impl() { slice_.shrink_to_fit(); }

		public:
		Mixed2D() = default;

		template<typename Points_rng>
		Mixed2D(Points_rng & points, const resolution_type res, const chs::flags::build::flags_t flags = {}) :
		        Mixed2D(points, dimensions_type(chs::n_tuple<Dim>(res)), flags)
		{}

		template<typename Points_rng>
		Mixed2D(Points_rng & points, const dimensions_type & res, const chs::flags::build::flags_t flags = {}) :
		        Super(points, res, flags), slice_(Box::mbb(points), res)
		{
			this->init(points, res, flags);
		}
	};
} // namespace chs