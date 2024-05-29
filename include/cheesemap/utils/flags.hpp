#pragma once

#include <cstddef>

namespace chs::flags::build
{
	enum
	{
		PARALLEL      = 1 << 0,
		REORDER       = 1 << 1,
		SHRINK_TO_FIT = 1 << 2,
	};

	using flags_t = std::size_t;
} // namespace chs::flags::build