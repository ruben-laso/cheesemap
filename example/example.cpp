#include <filesystem>
#include <iostream>

#include <cheesemap/cheesemap.hpp>

// abseil for flat_hash_map
#include <absl/container/flat_hash_map.h>

#include "handlers.hpp"

template<typename BuildMapFunction, typename Points>
void benchmark_map(const std::string & map_name, const BuildMapFunction && build_map, Points & points);

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

	const auto flags =
	        // chs::flags::build::REORDER | //
	        // chs::flags::build::PARALLEL | //
	        chs::flags::build::SHRINK_TO_FIT | //
	        0;

	benchmark_map("chs::Dense<2>", [&](auto & pts) { return chs::Dense<chs::Point, 2>(pts, 5.0, flags); }, points);
	benchmark_map("chs::Dense<3>", [&](auto & pts) { return chs::Dense<chs::Point, 3>(pts, 5.0, flags); }, points);

	benchmark_map(
	        "chs::Sparse<2, std::unordered_map>",
	        [&](auto & pts) { return chs::Sparse<chs::Point, 2, std::unordered_map>(pts, 5.0, flags); }, points);
	benchmark_map(
	        "chs::Sparse<3, std::unordered_map>",
	        [&](auto & pts) { return chs::Sparse<chs::Point, 3, std::unordered_map>(pts, 5.0, flags); }, points);

	benchmark_map(
	        "chs::Sparse<2, std::map>",
	        [&](auto & pts) { return chs::Sparse<chs::Point, 2, std::map>(pts, 5.0, flags); }, points);
	benchmark_map(
	        "chs::Sparse<3, std::map>",
	        [&](auto & pts) { return chs::Sparse<chs::Point, 3, std::map>(pts, 5.0, flags); }, points);

	benchmark_map(
	        "chs::Sparse<2, absl::flat_hash_map>",
	        [&](auto & pts) { return chs::Sparse<chs::Point, 2, absl::flat_hash_map>(pts, 5.0, flags); }, points);
	benchmark_map(
	        "chs::Sparse<3, absl::flat_hash_map>",
	        [&](auto & pts) { return chs::Sparse<chs::Point, 3, absl::flat_hash_map>(pts, 5.0, flags); }, points);

	benchmark_map("chs::Mixed<2>", [&](auto & pts) { return chs::Mixed<chs::Point, 2>(pts, 5.0, flags); }, points);
	benchmark_map("chs::Mixed<3>", [&](auto & pts) { return chs::Mixed<chs::Point, 3>(pts, 5.0, flags); }, points);

	return EXIT_SUCCESS;
}

template<typename BuildMapFunction, typename Points>
void benchmark_map(const std::string & map_name, const BuildMapFunction && build_map, Points & points)
{
	std::cout << "----------------------------------------\n";
	std::cout << "Benchmarking map: " << map_name << '\n';

	const auto start = std::chrono::high_resolution_clock::now();
	const auto map   = build_map(points);
	const auto end   = std::chrono::high_resolution_clock::now();
	std::cout << "chs:: map build time: " << std::chrono::duration<double>(end - start).count() << " s\n";

	const auto points_per_cell = map.points_per_cell();
	const auto non_empty_cells = ranges::count_if(points_per_cell, [](const auto & count) { return count > 0; });
	std::cout << "Number of cells: " << points_per_cell.size() << '\n';
	std::cout << "Number of non-empty cells: " << non_empty_cells << " ("
	          << (non_empty_cells * 100.0 / points_per_cell.size()) << "%)\n";

	// Count stored cells
	const auto stored_cells    = map.cells_stored();
	const auto no_stored_cells = ranges::count_if(stored_cells, [](const auto & exists) { return exists; });
	std::cout << "Stored cells: " << no_stored_cells << " (" << (no_stored_cells * 100.0 / points_per_cell.size())
	          << "%)\n";

	// Compute average points per cell
	const auto av_points_per_cell =
	        ranges::accumulate(points_per_cell, 0.0, std::plus<>()) / static_cast<double>(points_per_cell.size());
	std::cout << "Average points per cell: " << av_points_per_cell << '\n';
	// Average points per non-empty cell
	const auto av_points_per_non_empty_cell =
	        ranges::accumulate(points_per_cell, 0.0, std::plus<>()) / static_cast<double>(non_empty_cells);
	std::cout << "Average points per non-empty cell: " << av_points_per_non_empty_cell << '\n';

	benchmark_query(map, points);
	benchmark_knn(map, points);

	std::cout << "----------------------------------------\n";
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