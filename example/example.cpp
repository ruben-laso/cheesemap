#include <filesystem>
#include <iostream>

#include <cheesemap/cheesemap.hpp>

#include "handlers.hpp"

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
	auto map   = chs::Dense<chs::Point, 2>(points, 5.0);
	// auto map   = chs::Dense<chs::Point, 3>(points, 5.0);
	// auto map = chs::Dense2D<chs::Point>(points, 5.0);
	// auto map = chs::Dense3D<chs::Point>(points, 5.0);
	end        = std::chrono::high_resolution_clock::now();
	std::cout << "Dense2D map build time: "
	          << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms\n";

	std::size_t max_neighs = std::numeric_limits<std::size_t>::min();
	std::size_t min_neighs = std::numeric_limits<std::size_t>::max();
	double      ave_neighs = 0.0;

	std::size_t ns_map = 0;
	for (const auto & p : points | ranges::views::take(1'000))
	{
		chs::kernels::Sphere<3> search(p, 2.5);

		start                  = std::chrono::high_resolution_clock::now();
		const auto results_map = map.query(search);
		end                    = std::chrono::high_resolution_clock::now();

		ns_map += static_cast<std::size_t>(
		        std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

		max_neighs = std::max(max_neighs, results_map.size());
		min_neighs = std::min(min_neighs, results_map.size());
		ave_neighs += static_cast<double>(results_map.size()) / static_cast<double>(points.size());
	}

	std::cout << "Total search time: " << ns_map << "ns" << '\n';
	std::cout << "Average search time: " << ns_map / 10'000 << "ns" << '\n';
	std::cout << "Max neighbors: " << max_neighs << '\n';
	std::cout << "Min neighbors: " << min_neighs << '\n';
	std::cout << "Average neighbors: " << ave_neighs << '\n';

	return 0;
}
