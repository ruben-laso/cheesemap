#pragma once

#include <concepts>

#include "cheesemap/utils/Box.hpp"

namespace chs::concepts
{
	/**
	 * @brief A concept for a search kernel.
	 *
	 * @tparam Kernel_type The type of the kernel.
	 * @tparam Point_type The type of the point.	 *
	 */
	template<typename Kernel_type, typename Point_type>
	concept Kernel = requires(Kernel_type kernel, const Point_type & point, const chs::Box & box) {
		{
			kernel.is_inside(point)
		} -> std::convertible_to<bool>;
		{
			kernel.box()
		} -> std::convertible_to<const chs::Box &>;
	};
} // namespace chs::concepts
