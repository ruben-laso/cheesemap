#pragma once

#include <armadillo>
#include <utility>
#include <vector>

#include <range/v3/all.hpp>

template<typename... T>
inline bool are_the_same(const ranges::range auto & v1_, const ranges::range auto & v2_)
{
	if (v1_.size() != v2_.size())
	{
		std::cerr << "Vectors are not of the same length" << '\n';
		return false;
	}

	auto v1 = v1_ | ranges::to_vector;
	auto v2 = v2_ | ranges::to_vector;
	std::sort(std::begin(v1), std::end(v1));
	std::sort(std::begin(v2), std::end(v2));

	size_t mismatches = 0;

	for (const auto & [e1, e2] : ranges::views::zip(v1, v2))
	{
		if (e1 not_eq e2) { ++mismatches; }
	}

	if (mismatches > 0)
	{
		std::cerr << "Vectors have " << mismatches << " mismatches" << '\n'; //
	}

	return mismatches == 0;
}

template<typename T>
auto get_random(const T & min, const T & max)
{
	static std::random_device rd;
	static std::mt19937       gen(rd());

	if constexpr (std::is_floating_point_v<T>)
	{
		std::uniform_real_distribution<T> dis(min, max);
		return dis(gen);
	}
	else
	{
		std::uniform_int_distribution<T> dis(min, max);
		return dis(gen);
	}
}

auto random_bool()
{
	return std::cmp_equal(get_random(0, 1), 1);
}