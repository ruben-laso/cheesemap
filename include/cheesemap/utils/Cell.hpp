#pragma once

#include <vector>

namespace chs
{
	template<typename Point_type>
	using Cell = std::vector<Point_type *>;
} // namespace chs