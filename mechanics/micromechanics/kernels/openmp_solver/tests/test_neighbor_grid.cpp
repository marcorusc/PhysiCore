#include <algorithm>
#include <random>
#include <vector>

#include <gtest/gtest.h>

#include "neighbor_grid.h"

namespace physicore::mechanics::micromechanics::kernels::openmp_solver::tests {

namespace {
std::vector<real_t> random_positions(index_t agents_count, index_t dims, real_t lowest, real_t highest, unsigned seed)
{
	std::mt19937 generator(seed);
	std::uniform_real_distribution<real_t> coordinate(lowest, highest);

	std::vector<real_t> positions(agents_count * dims);
	std::ranges::generate(positions, [&] { return coordinate(generator); });
	return positions;
}

std::vector<index_t> candidates_of(const neighbor_grid& grid, const std::vector<real_t>& positions, index_t dims,
								   index_t i)
{
	std::vector<index_t> candidates;
	grid.visit_candidates(&positions[i * dims], [&candidates](index_t j) { candidates.push_back(j); });
	return candidates;
}

real_t squared_distance(const std::vector<real_t>& positions, index_t dims, index_t i, index_t j)
{
	real_t result = 0.0;
	for (index_t d = 0; d < dims; ++d)
	{
		const real_t difference = positions[i * dims + d] - positions[j * dims + d];
		result += difference * difference;
	}
	return result;
}

index_t bins_count(const neighbor_grid& grid) { return grid.shape()[0] * grid.shape()[1] * grid.shape()[2]; }
} // namespace

class NeighborGridTest : public ::testing::TestWithParam<index_t>
{};

INSTANTIATE_TEST_SUITE_P(Dimensions, NeighborGridTest, ::testing::Values(1, 2, 3));

TEST_P(NeighborGridTest, EveryPairWithinRangeIsFound)
{
	const index_t dims = GetParam();
	const index_t agents_count = 300;
	const real_t range = 1.0;
	// negative coordinates included; the box is several bins wide
	const auto positions = random_positions(agents_count, dims, -5.0, 7.0, 11);

	neighbor_grid grid;
	grid.build(positions.data(), agents_count, dims, range);

	index_t pairs_within_range = 0;
	for (index_t i = 0; i < agents_count; ++i)
	{
		auto candidates = candidates_of(grid, positions, dims, i);
		std::ranges::sort(candidates);

		// nobody is visited twice, and the sub-agent finds itself
		EXPECT_EQ(std::ranges::adjacent_find(candidates), candidates.end());
		EXPECT_TRUE(std::ranges::binary_search(candidates, i));

		for (index_t j = 0; j < agents_count; ++j)
		{
			if (i != j && squared_distance(positions, dims, i, j) < range * range)
			{
				++pairs_within_range;
				EXPECT_TRUE(std::ranges::binary_search(candidates, j)) << "pair " << i << ", " << j;
			}
		}
	}

	// the test would be empty if nobody had a neighbour
	EXPECT_GT(pairs_within_range, index_t { 100 });
}

TEST_P(NeighborGridTest, VisitsFarFewerCandidatesThanAllPairs)
{
	const index_t dims = GetParam();
	const index_t agents_count = 500;
	const auto positions = random_positions(agents_count, dims, 0.0, 40.0, 5);

	neighbor_grid grid;
	grid.build(positions.data(), agents_count, dims, 1.0);

	std::size_t visited = 0;
	for (index_t i = 0; i < agents_count; ++i)
	{
		visited += candidates_of(grid, positions, dims, i).size();
	}

	EXPECT_LT(visited, static_cast<std::size_t>(agents_count) * agents_count / 5);
}

TEST_P(NeighborGridTest, WidelySpreadSubAgentsDoNotBlowUpTheGrid)
{
	const index_t dims = GetParam();

	// two close sub-agents and two that are a billion interaction ranges away
	std::vector<real_t> positions(4 * dims, 0.0);
	positions[1 * dims] = 0.5;
	positions[2 * dims] = 1e9;
	positions[3 * dims] = -1e9;

	neighbor_grid grid;
	grid.build(positions.data(), 4, dims, 1.0);

	EXPECT_LE(bins_count(grid), index_t { 4 * 4 + 4096 });

	const auto candidates = candidates_of(grid, positions, dims, 0);
	EXPECT_NE(std::ranges::find(candidates, index_t { 1 }), candidates.end());
}

TEST_P(NeighborGridTest, CoincidentSubAgentsShareOneBin)
{
	const index_t dims = GetParam();
	const std::vector<real_t> positions(5 * dims, 3.25);

	neighbor_grid grid;
	grid.build(positions.data(), 5, dims, 1.0);

	EXPECT_EQ(bins_count(grid), index_t { 1 });
	EXPECT_EQ(candidates_of(grid, positions, dims, 2).size(), 5U);
}

TEST_P(NeighborGridTest, PairAtExactlyOneBinWidthIsFound)
{
	const index_t dims = GetParam();

	// sub-agents spaced by exactly the interaction range sit on bin boundaries
	std::vector<real_t> positions(6 * dims, 0.0);
	for (index_t i = 0; i < 6; ++i)
	{
		positions[i * dims] = static_cast<real_t>(i);
	}

	neighbor_grid grid;
	grid.build(positions.data(), 6, dims, 1.0);

	for (index_t i = 0; i + 1 < 6; ++i)
	{
		const auto candidates = candidates_of(grid, positions, dims, i);
		EXPECT_NE(std::ranges::find(candidates, i + 1), candidates.end()) << "sub-agent " << i;
	}
}

TEST(NeighborGridTest, EmptySystem)
{
	neighbor_grid grid;
	grid.build(nullptr, 0, 3, 1.0);

	EXPECT_EQ(bins_count(grid), index_t { 1 });

	const std::vector<real_t> probe = { 0.0, 0.0, 0.0 };
	EXPECT_TRUE(candidates_of(grid, probe, 3, 0).empty());
}

TEST(NeighborGridTest, PositionOutsideTheBuiltBoxUsesTheNearestBin)
{
	const std::vector<real_t> positions = { 0.0, 0.0, 10.0, 10.0 };

	neighbor_grid grid;
	grid.build(positions.data(), 2, 2, 1.0);

	const std::vector<real_t> far_below = { -50.0, -50.0 };
	const auto candidates = candidates_of(grid, far_below, 2, 0);

	ASSERT_EQ(candidates.size(), 1U);
	EXPECT_EQ(candidates[0], index_t { 0 });
}

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver::tests
