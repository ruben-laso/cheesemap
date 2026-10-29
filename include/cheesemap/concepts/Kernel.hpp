#pragma once

#include <concepts>

#include "cheesemap/utils/Box.hpp"

namespace chs::concepts
{
	/**
	 * @brief A search kernel that can test the map's point type.
	 *
	 * @tparam KernelType The type of the kernel.
	 * @tparam PointType The type of the point.
	 */
	template<typename KernelType, typename PointType>
	concept Kernel = requires(const KernelType & kernel, const PointType & point) {
		{ kernel.is_inside(point) } -> std::convertible_to<bool>;
		{ kernel.box() };
		{ KernelType::dim } -> std::convertible_to<std::size_t>;
	};
} // namespace chs::concepts
