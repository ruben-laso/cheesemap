#pragma once

#include <concepts>

namespace chs::concepts
{
	/**
	 * @brief A concept for a filter.
	 *
	 * @tparam Filter_type The type of the filter.
	 * @tparam Point_type The type of the point.
	 */
	template<typename Filter_type, typename Point_type>
	concept Filter = requires(Filter_type filter, const Point_type & point) {
		{
			filter(point)
		} -> std::convertible_to<bool>;
	};
} // namespace chs::concepts