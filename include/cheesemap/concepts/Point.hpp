#pragma once

#include <concepts>
#include <cstddef>
#include <type_traits>

namespace chs::concepts
{
	/**
	 * @brief A concept for a point in a multi-dimensional space.
	 *
	 * A point must be default constructible, copy constructible/assignable and indexable (read/write) with `operator[]`.
	 *
	 * @tparam P The type of the point (e.g., Eigen::Vector3d, arma::vec3, std::array<double, 3>).
	 * @tparam Dim The number of dimensions of the point.
	 */
	template<typename P, std::size_t Dim>
	concept Point = (Dim > 0) and std::default_initializable<P> and std::copy_constructible<P> and
	                std::is_copy_assignable_v<P> and requires(P & p, const P & cp, std::size_t i) {
		                { cp[i] } -> std::convertible_to<double>;
		                { p[i] } -> std::assignable_from<double>;
	                };
} // namespace chs::concepts
