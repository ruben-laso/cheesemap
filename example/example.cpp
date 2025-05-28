#include <filesystem>
#include <iostream>

#include <cheesemap/cheesemap.hpp>

#include "handlers.hpp"

template<typename Map, typename Points>
void benchmark_query(const Map & map, const Points & points);

template<typename Map, typename Points>
void benchmark_knn(const Map & map, const Points & points);

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

	start = std::chrono::high_resolution_clock::now();

	const auto flags = chs::flags::build::REORDER | chs::flags::build::PARALLEL | chs::flags::build::SHRINK_TO_FIT;

	static constexpr std::size_t Dims = 2;

	// auto map = chs::Dense<chs::Point, Dims>(points, 5.0, flags);
	// auto map = chs::Sparse<chs::Point, Dims>(points, 5.0, flags);
	auto map = chs::Mixed<chs::Point, Dims>(points, 5.0, flags);

	end = std::chrono::high_resolution_clock::now();
	std::cout << "chs:: map build time: "
	          << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms\n";

	const auto bytes = map.mem_footprint();
	const auto mb    = bytes / (1024.0 * 1024.0);
	std::cout << "Estimated mem. footprint: " << map.mem_footprint() << " Bytes (" << mb << "MB)" << '\n';

	const auto points_per_cell = map.points_per_cell();
	const auto non_empty_cells = ranges::count_if(points_per_cell, [](const auto & count) { return count > 0; });
	std::cout << "Number of cells: " << points_per_cell.size() << '\n';
	std::cout << "Number of non-empty cells: " << non_empty_cells << " (" << (non_empty_cells * 100.0 / points_per_cell.size()) << "%)\n";
	// Compute average points per cell
	const auto av_points_per_cell = ranges::fold_left(points_per_cell, 0.0, std::plus<>()) / static_cast<double>(points_per_cell.size());
	std::cout << "Average points per cell: " << av_points_per_cell << '\n';
	// Average points per non-empty cell
	const auto av_points_per_non_empty_cell = ranges::fold_left(points_per_cell, 0.0, std::plus<>()) / static_cast<double>(non_empty_cells);
	std::cout << "Average points per non-empty cell: " << av_points_per_non_empty_cell << '\n';

	benchmark_query(map, points);

	benchmark_knn(map, points);

	return EXIT_SUCCESS;
}

template<typename Map, typename Points>
void benchmark_query(const Map & map, const Points & points)
{
	std::size_t max_neighs = std::numeric_limits<std::size_t>::min();
	std::size_t min_neighs = std::numeric_limits<std::size_t>::max();
	double      ave_neighs = 0.0;

	static constexpr std::size_t NUM_SEARCHES = 100'000;

	std::size_t ns_map = 0;
	for (const auto & p : points | ranges::views::take(NUM_SEARCHES))
	{
		chs::kernels::Sphere<3> search(p, 2.5);

		const auto start       = std::chrono::high_resolution_clock::now();
		const auto results_map = map.query(search);
		const auto end         = std::chrono::high_resolution_clock::now();

		ns_map += static_cast<std::size_t>(
		        std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

		max_neighs = std::max(max_neighs, results_map.size());
		min_neighs = std::min(min_neighs, results_map.size());
		ave_neighs += static_cast<double>(results_map.size()) / static_cast<double>(points.size());
	}

	std::cout << "#Queries: " << NUM_SEARCHES << '\n';
	std::cout << "Total query search time: " << ns_map << "ns" << '\n';
	std::cout << "Average query search time: " << ns_map / NUM_SEARCHES << "ns" << '\n';
	std::cout << "Max neighbors: " << max_neighs << '\n';
	std::cout << "Min neighbors: " << min_neighs << '\n';
	std::cout << "Average neighbors: " << ave_neighs << '\n';
}

template<typename Map, typename Points>
void benchmark_knn(const Map & map, const Points & points)
{
	static constexpr std::size_t NUM_SEARCHES = 100'000;

	static constexpr std::size_t NUM_KNN = 10;

	std::size_t ns_map = 0;
	for (const auto & p : points | ranges::views::take(NUM_SEARCHES))
	{
		const auto start       = std::chrono::high_resolution_clock::now();
		const auto results_map = map.knn(NUM_KNN, p);
		const auto end         = std::chrono::high_resolution_clock::now();

		ns_map += static_cast<std::size_t>(
		        std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
	}

	std::cout << "#knn searches: " << NUM_SEARCHES << '\n';
	std::cout << "Total knn search time: " << ns_map << "ns" << '\n';
	std::cout << "Average knn search time: " << ns_map / NUM_SEARCHES << "ns" << '\n';
}