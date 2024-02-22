#pragma once

#include <memory>
#include <vector>

#include <range/v3/all.hpp>

#include "Box.hpp"

namespace chs
{
	template<typename Point_type>
	using Cell = std::vector<Point_type *>;
} // namespace chs