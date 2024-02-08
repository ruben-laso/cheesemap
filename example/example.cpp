#include <filesystem>
#include <iostream>

#include <cheesemap/cheesemap.hpp>

#include "handlers.hpp"
#include "Spherical_search.hpp"

template<typename... T>
inline bool are_the_same(const std::vector<T...> & v1_, const std::vector<T...> & v2_)
{
	if (v1_.size() != v2_.size())
	{
		std::cerr << "Vectors are not of the same length" << '\n';
		return false;
	}

	auto v1 = v1_;
	auto v2 = v2_;
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

auto ground_truth(std::vector<chs::Point> & points, const chs::Spherical_search<3> & search)
        -> std::vector<chs::Point *>
{
	return ranges::views::filter(points, [&](const auto & neigh) { return search.is_inside(neigh); }) |
	       ranges::views::transform([](auto & p) { return &p; }) | ranges::to<std::vector>;
}

auto main(const int argc, const char * const argv[]) -> int
{
	if (argc < 2)
	{
		std::cerr << "Usage: " << argv[0] << " <path_to_map>" << std::endl;
		return 1;
	}

	auto start = std::chrono::high_resolution_clock::now();
	auto end   = std::chrono::high_resolution_clock::now();
	std::cout << "Clock resolution: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()
	          << "ns\n";

	std::filesystem::path path = argv[1];

	start       = std::chrono::high_resolution_clock::now();
	auto points = readPointCloud(path);
	end         = std::chrono::high_resolution_clock::now();
	std::cout << "Reading input file time: "
	          << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms\n";

	start      = std::chrono::high_resolution_clock::now();
	auto map2d = chs::dense2d<chs::Point>(points, 1.0);
	end        = std::chrono::high_resolution_clock::now();
	std::cout << "dense2d map build time: "
	          << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms\n";

	std::size_t max_neighs = std::numeric_limits<std::size_t>::min();
	std::size_t min_neighs = std::numeric_limits<std::size_t>::max();
	double      ave_neighs = 0.0;

	std::size_t ns_map          = 0;
	std::size_t ns_ground_truth = 0;
	for (const auto & p : points | ranges::views::take(1'000))
	{
		chs::Spherical_search<3> search(p, 1.0);

		start                  = std::chrono::high_resolution_clock::now();
		const auto results_map = map2d.query(search);
		end                    = std::chrono::high_resolution_clock::now();

		ns_map += static_cast<std::size_t>(
		        std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

		start                    = std::chrono::high_resolution_clock::now();
		const auto results_truth = ground_truth(points, search);
		end                      = std::chrono::high_resolution_clock::now();

		ns_ground_truth += static_cast<std::size_t>(
		        std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

		if (not are_the_same(results_map, results_truth))
		{
			std::cerr << "Results of map and ground truth are not the same" << '\n';
			// Repeat computation to debug
			std::ignore = map2d.query(search);
			return 1;
		}
		else
		{
			// std::cout << "Results of map and ground truth are the same" << '\n';
		}

		max_neighs = std::max(max_neighs, results_map.size());
		min_neighs = std::min(min_neighs, results_map.size());
		ave_neighs += static_cast<double>(results_map.size()) / static_cast<double>(points.size());
	}

	std::cout << "Total search time (map / ground truth): " << ns_map << "ns / " << ns_ground_truth << "ns\n";
	std::cout << "Average search time (map / ground truth): " << ns_map / 10'000 << "ns / "
	          << ns_ground_truth / 10'000 << "ns\n";
	std::cout << "Max neighbors: " << max_neighs << '\n';
	std::cout << "Min neighbors: " << min_neighs << '\n';
	std::cout << "Average neighbors: " << ave_neighs << '\n';

	return 0;
}
