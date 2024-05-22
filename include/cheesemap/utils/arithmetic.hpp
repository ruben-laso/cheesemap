#pragma once

#include <algorithm>
#include <utility>

#include <cmath>

namespace chs
{
	template<std::size_t... Is>
	[[nodiscard]] inline auto sq_distance(const auto & p1, const auto & p2, std::index_sequence<Is...>) -> double
	{
		using T = std::common_type_t<decltype(p1[Is])...>;

		T dist = 0;

		// This is a fold expression: dist^2 = \sum_i (p1[i] - p2[i])^2
		((dist += (p1[Is] - p2[Is]) * (p1[Is] - p2[Is])), ...);

		return dist;
	}

	template<std::size_t Dim = 3>
	[[nodiscard]] inline auto sq_distance(const auto & p1, const auto & p2) -> double
	{
		return sq_distance(p1, p2, std::make_index_sequence<Dim>{});
	}

	template<std::size_t Dim = 3>
	[[nodiscard]] inline auto distance(const auto & p1, const auto & p2) -> double
	{
		return std::sqrt(sq_distance(p1, p2, std::make_index_sequence<Dim>{}));
	}

	template<std::size_t... Is>
	inline void clamp(auto & val, const auto & min, const auto & max, std::index_sequence<Is...>)
	{
		((val[Is] = std::clamp(val[Is], min[Is], max[Is])), ...);
	}

	template<std::size_t Dim>
	inline void clamp(auto & vals, const auto & mins, const auto & maxs)
	{
		chs::clamp(vals, mins, maxs, std::make_index_sequence<Dim>{});
	}

	[[nodiscard]] inline auto radius_for_density(const auto & curr_pts, const auto & curr_r, const auto & trgt_pts)
	{
		const auto trgt_radius =
		        curr_r * std::cbrt(static_cast<double>(trgt_pts) / static_cast<double>(curr_pts));

		return trgt_radius;
	}

	template<std::integral T, std::integral U, std::integral V>
	[[nodiscard]] inline auto within_closed_bounds(const T & val, const U & min, const V & max) -> bool
	{
		return std::cmp_less_equal(min, val) and std::cmp_less_equal(val, max);
	}

	[[nodiscard]] inline auto within_closed_bounds(const auto & val_min_max) -> bool
	{
		const auto & [val, min, max] = val_min_max;
		return within_closed_bounds(val, min, max);
	}

	template<std::size_t... Is>
	[[nodiscard]] inline auto within_closed_bounds(const auto & vals, const auto & mins, const auto & maxs,
	                                               std::index_sequence<Is...>)
	{
		return (within_closed_bounds(vals[Is], mins[Is], maxs[Is]) and ...);
	}

	template<std::size_t Dim = 3>
	[[nodiscard]] inline auto within_closed_bounds(const auto & vals, const auto & mins, const auto & maxs)
	{
		return within_closed_bounds(vals, mins, maxs, std::make_index_sequence<Dim>{});
	}

	template<std::size_t... Is>
	[[nodiscard]] inline auto all_visited(const auto & mins, const auto & maxs, const auto & sizes,
	                                      std::index_sequence<Is...>)
	{
		return (std::cmp_greater_equal(maxs[Is] - mins[Is], sizes[Is] - 1) and ...);
	}

	template<std::size_t Dim = 3>
	[[nodiscard]] inline auto all_visited(const auto & mins, const auto & maxs, const auto & sizes)
	{
		return all_visited(mins, maxs, sizes, std::make_index_sequence<Dim>{});
	}

	template<std::size_t... Is>
	[[nodiscard]] inline auto all_equal(const auto & as, const auto & bs, std::index_sequence<Is...>)
	{
		return ((as[Is] == bs[Is]) and ...);
	}

	template<std::size_t Dim>
	[[nodiscard]] inline auto all_equal(const auto & as, const auto & bs)
	{
		return all_equal(as, bs, std::make_index_sequence<Dim>{});
	}
} // namespace chs