#pragma once

#include <array>
#include <utility>

#include <range/v3/view/cartesian_product.hpp>
#include <range/v3/view/indices.hpp>
#include <range/v3/view/transform.hpp>

namespace chs
{
	template<typename RangeA, typename RangeB, std::size_t... Is>
	constexpr auto cartesian_product_helper(RangeA && mins, RangeB && maxs, std::index_sequence<Is...>)
	{
		return ranges::views::cartesian_product(
		        ranges::views::closed_indices(std::get<Is>(mins), std::get<Is>(maxs))...);
	}

	template<std::size_t N, typename RangeA, typename RangeB>
	constexpr auto cartesian(RangeA && mins, RangeB && maxs)
	{
		return cartesian_product_helper(mins, maxs, std::make_index_sequence<N>{});
	}
	template<typename RangeA, typename RangeB, std::size_t... Is>
	constexpr auto cartesian_product_helper_open_bound(RangeA && mins, RangeB && maxs, std::index_sequence<Is...>)
	{
		return ranges::views::cartesian_product(
		        ranges::views::indices(std::get<Is>(mins), std::get<Is>(maxs))...);
	}

	template<std::size_t N, typename RangeA, typename RangeB>
	constexpr auto cartesian_open_bound(RangeA && mins, RangeB && maxs)
	{
		return cartesian_product_helper_open_bound(mins, maxs, std::make_index_sequence<N>{});
	}

	inline auto tuple_to_array(const auto & tuple)
	{
		return std::apply([](auto... args) { return std::array{ args... }; }, tuple);
	}

	template<std::size_t N, typename Range>
	constexpr auto cartesian_as_array(Range && mins, Range && maxs)
	{
		return cartesian<N>(mins, maxs) |
		       ranges::views::transform([](const auto & tuple) { return tuple_to_array(tuple); });
	}
} // namespace chs