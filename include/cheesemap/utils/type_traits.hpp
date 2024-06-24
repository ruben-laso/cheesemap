#pragma once

#include <tuple>

namespace chs::type_traits
{
	// From https://stackoverflow.com/a/33511913
	template<typename T, unsigned N, typename... REST>
	struct generate_tuple_type
	{
		typedef typename generate_tuple_type<T, N-1, T, REST...>::type type;
	};

	template<typename T, typename... REST>
	struct generate_tuple_type<T, 0, REST...>
	{
		typedef std::tuple<REST...> type;
	};

	template<typename T, unsigned N>
	using tuple = typename generate_tuple_type<T, N>::type;
} // namespace chs::type_traits