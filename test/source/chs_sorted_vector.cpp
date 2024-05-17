#include <cheesemap/utils/sorted_vector.hpp>

#include <gtest/gtest.h>

#include <range/v3/view/indices.hpp>

#include <random>


template<std::integral int_t = int>
auto rand_int() -> int_t
{
	static std::random_device                   rd;
	static std::mt19937                         gen(rd());
	static std::uniform_int_distribution<int_t> dis(std::numeric_limits<int_t>::min(),
	                                                std::numeric_limits<int_t>::max());
	return dis(gen);
}

TEST(chs, sorted_vector)
{
	static constexpr std::size_t NUM_ELEMS = 10;

	chs::sorted_vector<int> sorted_vector;

	EXPECT_TRUE(sorted_vector.empty());
	EXPECT_EQ(sorted_vector.size(), 0);

	// Insert some elements
	for ([[maybe_unused]] const auto _ : ranges::views::indices(NUM_ELEMS))
	{
		sorted_vector.insert(rand_int());
	}

	EXPECT_FALSE(sorted_vector.empty());
	EXPECT_EQ(sorted_vector.size(), NUM_ELEMS);
	EXPECT_TRUE(std::is_sorted(sorted_vector.begin(), sorted_vector.end()));
}

TEST(chs, sorted_vector_max_size)
{
	static constexpr std::size_t NUM_ELEMS = 10;

	chs::sorted_vector<int> sorted_vector(NUM_ELEMS);

	EXPECT_TRUE(sorted_vector.empty());
	EXPECT_EQ(sorted_vector.size(), 0);

	// Insert some elements
	for ([[maybe_unused]] const auto _ : ranges::views::indices(2 * NUM_ELEMS))
	{
		sorted_vector.insert(rand_int());
	}

	EXPECT_FALSE(sorted_vector.empty());
	EXPECT_EQ(sorted_vector.size(), NUM_ELEMS);
	EXPECT_TRUE(std::is_sorted(sorted_vector.begin(), sorted_vector.end()));
}

TEST(chs, sorted_vector_compare)
{
	static constexpr std::size_t NUM_ELEMS = 10;

	chs::sorted_vector<int, std::greater<int>> sorted_vector;

	EXPECT_TRUE(sorted_vector.empty());
	EXPECT_EQ(sorted_vector.size(), 0);

	// Insert some elements
	for ([[maybe_unused]] const auto _ : ranges::views::indices(NUM_ELEMS))
	{
		sorted_vector.insert(rand_int());
	}

	EXPECT_FALSE(sorted_vector.empty());
	EXPECT_EQ(sorted_vector.size(), NUM_ELEMS);
	EXPECT_TRUE(std::is_sorted(sorted_vector.begin(), sorted_vector.end(), std::greater<int>()));
}

TEST(chs, sorted_vector_compare_max_size)
{
	static constexpr std::size_t NUM_ELEMS = 10;

	chs::sorted_vector<int, std::greater<int>> sorted_vector(NUM_ELEMS);

	EXPECT_TRUE(sorted_vector.empty());
	EXPECT_EQ(sorted_vector.size(), 0);

	// Insert some elements
	for ([[maybe_unused]] const auto _ : ranges::views::indices(2 * NUM_ELEMS))
	{
		sorted_vector.insert(rand_int());
	}

	EXPECT_FALSE(sorted_vector.empty());
	EXPECT_EQ(sorted_vector.size(), NUM_ELEMS);
	EXPECT_TRUE(std::is_sorted(sorted_vector.begin(), sorted_vector.end(), std::greater<int>()));
}

TEST(chs, sorted_vector_fuzzy)
{
	static constexpr std::size_t NUM_TESTS = 100'000;
	static constexpr std::size_t MAX_ELEMS = 50;

	for ([[maybe_unused]] const auto _ : ranges::views::indices(NUM_TESTS))
	{
		const auto num_elements = rand_int<std::size_t>() % MAX_ELEMS;

		chs::sorted_vector<int> sorted_vector(num_elements);

		EXPECT_TRUE(sorted_vector.empty());

		for ([[maybe_unused]] const auto i : ranges::views::indices(2 * num_elements))
		{
			sorted_vector.insert(rand_int());
		}

		EXPECT_EQ(sorted_vector.size(), num_elements);
		EXPECT_TRUE(std::is_sorted(sorted_vector.begin(), sorted_vector.end()));
	}
}

auto main() -> int
{
	::testing::InitGoogleTest();
	return RUN_ALL_TESTS();
}