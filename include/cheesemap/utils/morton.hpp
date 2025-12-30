#pragma once

#include "cheesemap/utils/Box.hpp"
#include "cheesemap/utils/Point.hpp"

namespace chs::morton
{
	namespace details
	{
		static constexpr auto BIT_DEPTH = 21; // 3 * 21 = 63 bits total
		static constexpr auto MAX_INT   = (1ULL << BIT_DEPTH) - 1;

		namespace Dim2
		{
			[[nodiscard]] inline auto expand_bits(std::uint64_t x) -> std::uint64_t
			{
				x &= 0xFFFFFFFF;
				x = (x | (x << 16)) & 0x0000FFFF0000FFFF;
				x = (x | (x << 8)) & 0x00FF00FF00FF00FF;
				x = (x | (x << 4)) & 0x0F0F0F0F0F0F0F0F;
				x = (x | (x << 2)) & 0x3333333333333333;
				x = (x | (x << 1)) & 0x5555555555555555;
				return x;
			}

			[[nodiscard]] inline auto morton(const std::uint64_t x, const std::uint64_t y) -> std::uint64_t
			{
				const auto xx = expand_bits(x);
				const auto yy = expand_bits(y);

				return (xx << 1) | yy;
			}

			[[nodiscard]] inline auto morton(const Point & p, const Box & bbox) -> std::uint64_t
			{
				const auto scale_x = (p[0] - bbox.min()[0]) / (bbox.max()[0] - bbox.min()[0]);
				const auto scale_y = (p[1] - bbox.min()[1]) / (bbox.max()[1] - bbox.min()[1]);

				const auto x = static_cast<std::uint64_t>(scale_x * MAX_INT);
				const auto y = static_cast<std::uint64_t>(scale_y * MAX_INT);

				return morton(x, y);
			}
		} // namespace Dim2

		namespace Dim3
		{
			static constexpr auto BIT_DEPTH = 21; // 3 * 21 = 63 bits total
			static constexpr auto MAX_INT   = (1ULL << BIT_DEPTH) - 1;

			[[nodiscard]] inline auto expand_bits(std::uint64_t x) -> std::uint64_t
			{
				x &= 0x1FFFFF;
				x = (x | (x << 32)) & 0x1F00000000FFFF;
				x = (x | (x << 16)) & 0x1F0000FF0000FF;
				x = (x | (x << 8)) & 0x100F00F00F00F00F;
				x = (x | (x << 4)) & 0x10C30C30C30C30C3;
				x = (x | (x << 2)) & 0x1249249249249249;
				return x;
			}

			[[nodiscard]] inline auto morton(const std::uint64_t x, const std::uint64_t y,
			                                 const std::uint64_t z) -> std::uint64_t
			{
				const auto xx = expand_bits(x);
				const auto yy = expand_bits(y);
				const auto zz = expand_bits(z);

				return (xx << 2) | (yy << 1) | zz;
			}

			[[nodiscard]] inline auto morton(const Point & p, const Box & bbox) -> std::uint64_t
			{
				const auto scale_x = (p[0] - bbox.min()[0]) / (bbox.max()[0] - bbox.min()[0]);
				const auto scale_y = (p[1] - bbox.min()[1]) / (bbox.max()[1] - bbox.min()[1]);
				const auto scale_z = (p[2] - bbox.min()[2]) / (bbox.max()[2] - bbox.min()[2]);

				const auto x = static_cast<std::uint64_t>(scale_x * MAX_INT);
				const auto y = static_cast<std::uint64_t>(scale_y * MAX_INT);
				const auto z = static_cast<std::uint64_t>(scale_z * MAX_INT);

				return morton(x, y, z);
			}
		} // namespace Dim3
	} // namespace details

	template<std::size_t Dim>
	[[nodiscard]] inline auto morton(const Point & p, const Box & bbox) -> std::uint64_t
	{
		static_assert(Dim == 2 or Dim == 3, "Only 2D and 3D Morton codes are supported.");
		if constexpr (Dim == 2) { return details::Dim2::morton(p, bbox); }
		else if (Dim == 3) { return details::Dim3::morton(p, bbox); }
	}
} // namespace chs::morton