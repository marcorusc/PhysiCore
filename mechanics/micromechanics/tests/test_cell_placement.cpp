#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "micromechanics/cell_placement.h"
#include "micromechanics/environment.h"

namespace physicore::mechanics::micromechanics::tests {

namespace {
real_t distance(const std::vector<real_t>& positions, index_t dims, index_t i, index_t j)
{
	real_t squared = 0.0;
	for (index_t d = 0; d < dims; ++d)
	{
		squared += std::pow(positions[i * dims + d] - positions[j * dims + d], 2);
	}
	return std::sqrt(squared);
}

real_t distance_from_origin(const std::vector<real_t>& positions, index_t dims, index_t i)
{
	real_t squared = 0.0;
	for (index_t d = 0; d < dims; ++d)
	{
		squared += std::pow(positions[i * dims + d], 2);
	}
	return std::sqrt(squared);
}
} // namespace

class CellPlacementTest : public ::testing::TestWithParam<index_t>
{};

INSTANTIATE_TEST_SUITE_P(Dimensions, CellPlacementTest, ::testing::Values(1, 2, 3));

TEST_P(CellPlacementTest, NoTwoElementsAreCloserThanTheSpacingAndNoneIsIsolated)
{
	const index_t dims = GetParam();
	const real_t spacing = 2.85;

	for (const index_t count : { 2, 7, 13, 50, 300 })
	{
		const auto positions = packed_ball_positions(dims, count, spacing);
		ASSERT_EQ(positions.size(), count * dims);

		for (index_t i = 0; i < count; ++i)
		{
			real_t nearest = std::numeric_limits<real_t>::max();
			for (index_t j = 0; j < count; ++j)
			{
				if (i != j)
				{
					nearest = std::min(nearest, distance(positions, dims, i, j));
				}
			}
			// exactly one spacing: closer would overlap, further would leave the element on its own
			EXPECT_NEAR(nearest, spacing, 1e-12) << count << " elements, element " << i;
		}
	}
}

TEST_P(CellPlacementTest, ElementsAreOrderedByDistanceStartingAtTheOrigin)
{
	const index_t dims = GetParam();
	const auto positions = packed_ball_positions(dims, 200, 1.0);

	EXPECT_DOUBLE_EQ(distance_from_origin(positions, dims, 0), 0.0);
	for (index_t i = 1; i < 200; ++i)
	{
		EXPECT_GE(distance_from_origin(positions, dims, i) + 1e-12, distance_from_origin(positions, dims, i - 1));
	}
}

TEST_P(CellPlacementTest, SmallerCellIsTheBeginningOfALargerOne)
{
	const index_t dims = GetParam();
	const auto small_cell = packed_ball_positions(dims, 40, 1.5);
	const auto large_cell = packed_ball_positions(dims, 500, 1.5);

	ASSERT_EQ(small_cell.size(), 40 * dims);
	EXPECT_TRUE(std::equal(small_cell.begin(), small_cell.end(), large_cell.begin()));
}

TEST_P(CellPlacementTest, CentreOfMassIsAtTheOriginForCompleteShells)
{
	const index_t dims = GetParam();
	// the origin and its complete first shell of nearest neighbours
	const std::array<index_t, 3> first_shell = { 2, 6, 12 };
	const index_t count = 1 + first_shell.at(dims - 1);

	const auto positions = packed_ball_positions(dims, count, 1.0);

	for (index_t d = 0; d < dims; ++d)
	{
		real_t sum = 0.0;
		for (index_t i = 0; i < count; ++i)
		{
			sum += positions[i * dims + d];
		}
		EXPECT_NEAR(sum, 0.0, 1e-12);
	}
	// every neighbour in the shell is exactly one spacing from the origin
	for (index_t i = 1; i < count; ++i)
	{
		EXPECT_NEAR(distance_from_origin(positions, dims, i), 1.0, 1e-12);
	}
}

TEST_P(CellPlacementTest, EmptyAndSingleElementCells)
{
	const index_t dims = GetParam();

	EXPECT_TRUE(packed_ball_positions(dims, 0, 1.0).empty());
	EXPECT_EQ(packed_ball_positions(dims, 1, 1.0), std::vector<real_t>(dims, 0.0));
}

TEST(CellPlacementTest, ThreeDimensionalCellFillsTheBallOfTheScalingRelation)
{
	// N elements at spacing r_eq = 2 R (p / N)^(1/3) with p the close-packing density fill a ball of
	// radius R (Sandersius and Newman 2008, Eq. 8). The lattice used here is close packed.
	const index_t count = 2000;
	const real_t radius = 10.0;
	const real_t packing_density = std::numbers::pi / (3.0 * std::numbers::sqrt2);
	const real_t spacing = 2.0 * radius * std::cbrt(packing_density / static_cast<real_t>(count));

	const auto positions = packed_ball_positions(3, count, spacing);

	real_t furthest = 0.0;
	real_t squared_gyration_radius = 0.0;
	for (index_t i = 0; i < count; ++i)
	{
		const real_t r = distance_from_origin(positions, 3, i);
		furthest = std::max(furthest, r);
		squared_gyration_radius += r * r / static_cast<real_t>(count);
	}

	// the outermost element centres lie on the surface of that ball, give or take one spacing
	EXPECT_NEAR(furthest, radius, spacing);
	// and are spread like a uniform ball of that radius
	EXPECT_NEAR(std::sqrt(5.0 / 3.0 * squared_gyration_radius), radius, 0.03 * radius);
}

TEST(CellPlacementTest, RejectsInvalidArguments)
{
	EXPECT_THROW(packed_ball_positions(0, 10, 1.0), std::invalid_argument);
	EXPECT_THROW(packed_ball_positions(4, 10, 1.0), std::invalid_argument);
	EXPECT_THROW(packed_ball_positions(3, 10, 0.0), std::invalid_argument);
	EXPECT_THROW(packed_ball_positions(3, 10, -1.0), std::invalid_argument);
}

TEST(CellPlacementTest, EnvironmentCreatesAPackedCellAroundACentre)
{
	environment env({ 3, { -100, -100, -100 }, { 100, 100, 100 }, { 20, 20, 20 } }, 0.1);
	env.params.equilibrium_distance = 2.0;
	const std::array<real_t, 3> centre = { 5.0, -3.0, 1.0 };

	env.create_cell();
	auto* cell = env.create_cell(centre, 13);

	ASSERT_EQ(env.cells->size(), 2U);
	ASSERT_EQ(env.agents->size(), 13U);
	EXPECT_EQ(cell->agent_count(), index_t { 13 });
	EXPECT_DOUBLE_EQ(cell->agent_radius(), 1.0);

	// the first sub-agent sits on the centre, the other twelve one equilibrium distance away
	for (index_t i = 0; i < 13; ++i)
	{
		auto* element = env.agents->get_agent_at(i);
		EXPECT_EQ(element->cell_index(), index_t { 1 });

		real_t squared = 0.0;
		for (std::size_t d = 0; d < 3; ++d)
		{
			squared += std::pow(element->position()[d] - centre.at(d), 2);
		}
		EXPECT_NEAR(std::sqrt(squared), i == 0 ? 0.0 : 2.0, 1e-12);
	}

	// a complete shell is symmetric, so the cell is at the requested centre
	for (std::size_t d = 0; d < 3; ++d)
	{
		EXPECT_NEAR(cell->position()[d], centre.at(d), 1e-12);
	}
}

TEST(CellPlacementTest, EnvironmentChecksTheCentre)
{
	environment env({ 2, { -100, -100, -100 }, { 100, 100, 100 }, { 20, 20, 20 } }, 0.1);
	const std::array<real_t, 3> too_many_values = { 0.0, 0.0, 0.0 };
	const std::array<real_t, 2> centre = { 1.0, 2.0 };

	EXPECT_THROW(env.create_cell(too_many_values, 5), std::invalid_argument);
	EXPECT_EQ(env.cells->size(), 0U);

	// in two dimensions the origin has six neighbours
	auto* cell = env.create_cell(centre, 7);
	EXPECT_EQ(cell->agent_count(), index_t { 7 });
	EXPECT_NEAR(cell->position()[0], 1.0, 1e-12);
	EXPECT_NEAR(cell->position()[1], 2.0, 1e-12);
}

} // namespace physicore::mechanics::micromechanics::tests
