#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "micromechanics/environment.h"

namespace physicore::mechanics::micromechanics::tests {

namespace {
environment make_environment(index_t dims)
{
	return { { dims, { -100, -100, -100 }, { 100, 100, 100 }, { 20, 20, 20 } }, 0.1 };
}

std::vector<real_t> positions_of(environment& env)
{
	std::vector<real_t> positions;
	for (index_t i = 0; i < env.agents->size(); ++i)
	{
		const auto position = env.agents->get_agent_at(i)->position();
		positions.insert(positions.end(), position.begin(), position.end());
	}
	return positions;
}

std::vector<index_t> cell_indices_of(environment& env)
{
	std::vector<index_t> indices;
	for (index_t i = 0; i < env.agents->size(); ++i)
	{
		indices.push_back(env.agents->get_agent_at(i)->cell_index());
	}
	return indices;
}

real_t distance(environment& env, index_t i, index_t j)
{
	const auto a = env.agents->get_agent_at(i)->position();
	const auto b = env.agents->get_agent_at(j)->position();

	real_t squared = 0.0;
	for (std::size_t d = 0; d < a.size(); ++d)
	{
		squared += std::pow(a[d] - b[d], 2);
	}
	return std::sqrt(squared);
}

// Distance from sub-agent i to the nearest other sub-agent
real_t nearest_neighbour(environment& env, index_t i)
{
	real_t nearest = std::numeric_limits<real_t>::max();
	for (index_t j = 0; j < env.agents->size(); ++j)
	{
		if (i != j)
		{
			nearest = std::min(nearest, distance(env, i, j));
		}
	}
	return nearest;
}

real_t projection(environment& env, index_t i, const std::vector<real_t>& normal)
{
	const auto position = env.agents->get_agent_at(i)->position();

	real_t sum = 0.0;
	for (std::size_t d = 0; d < position.size(); ++d)
	{
		sum += position[d] * normal[d];
	}
	return sum;
}
} // namespace

// ---------------------------------------------------------------------------------------------- growth

class GrowthTest : public ::testing::TestWithParam<index_t>
{};

INSTANTIATE_TEST_SUITE_P(Dimensions, GrowthTest, ::testing::Values(2, 3));

TEST_P(GrowthTest, NewSubAgentJoinsTheCellOneEquilibriumDistanceFromAMember)
{
	const index_t dims = GetParam();
	environment env = make_environment(dims);
	env.cell_types[0].equilibrium_distance = 2.0;
	const std::vector<real_t> centre(dims, 5.0);
	auto* cell = env.create_cell(centre, 13);

	auto* added = env.grow_cell(0);

	EXPECT_EQ(env.agents->size(), 14U);
	EXPECT_EQ(cell->agent_count(), index_t { 14 });
	EXPECT_EQ(added->cell_index(), index_t { 0 });
	EXPECT_EQ(env.agents->get_agent_at(13), added);
	for (const real_t v : added->velocity())
	{
		EXPECT_DOUBLE_EQ(v, 0.0);
	}

	// it touches at least one member of the cell and overlaps with none
	EXPECT_NEAR(nearest_neighbour(env, 13), 2.0, 1e-9);
}

TEST_P(GrowthTest, GrowingACellNeverCreatesStrongOverlaps)
{
	const index_t dims = GetParam();
	environment env = make_environment(dims);
	const std::vector<real_t> centre(dims, 0.0);
	env.create_cell(centre, 20);

	// without any relaxation in between, which is the hardest case
	real_t closest = std::numeric_limits<real_t>::max();
	for (index_t n = 0; n < 60; ++n)
	{
		env.grow_cell(0);
		closest = std::min(closest, nearest_neighbour(env, env.agents->size() - 1));
	}

	EXPECT_EQ(env.agents->size(), 80U);
	// no new sub-agent was placed closer than 0.7 equilibrium distances to an existing one
	EXPECT_GT(closest, 0.7);
}

TEST(GrowthTest, FirstSubAgentOfAnEmptyCellGoesToTheCell)
{
	environment env = make_environment(3);
	auto* cell = env.create_cell();
	cell->position()[0] = 4.0;
	cell->position()[2] = -1.0;

	auto* added = env.grow_cell(0);

	EXPECT_EQ(cell->agent_count(), index_t { 1 });
	EXPECT_DOUBLE_EQ(added->position()[0], 4.0);
	EXPECT_DOUBLE_EQ(added->position()[1], 0.0);
	EXPECT_DOUBLE_EQ(added->position()[2], -1.0);
}

TEST(GrowthTest, GrowthAvoidsNeighbouringCells)
{
	// Two cells in contact: the room around a candidate counts sub-agents of every cell.
	environment env = make_environment(3);
	const std::array<real_t, 3> left = { 0.0, 0.0, 0.0 };
	const std::array<real_t, 3> right = { 3.0, 0.0, 0.0 };
	env.create_cell(left, 13);
	env.create_cell(right, 13);

	for (index_t n = 0; n < 20; ++n)
	{
		env.grow_cell(0);
		EXPECT_GT(nearest_neighbour(env, env.agents->size() - 1), 0.7);
	}
	EXPECT_EQ(env.cells->get_agent_at(0)->agent_count(), index_t { 33 });
	EXPECT_EQ(env.cells->get_agent_at(1)->agent_count(), index_t { 13 });
}

TEST(GrowthTest, SameSeedGrowsTheSameCellAndAnotherSeedADifferentOne)
{
	const auto grow = [](std::uint64_t seed) {
		environment env = make_environment(3);
		env.params.random_seed = seed;
		const std::array<real_t, 3> centre = { 0.0, 0.0, 0.0 };
		env.create_cell(centre, 13);
		for (index_t n = 0; n < 30; ++n)
		{
			env.grow_cell(0);
		}
		return positions_of(env);
	};

	const auto first = grow(3);
	EXPECT_EQ(grow(3), first);
	EXPECT_NE(grow(4), first);
}

TEST(GrowthTest, RejectsAnUnknownCell)
{
	environment env = make_environment(3);
	EXPECT_THROW(env.grow_cell(0), std::out_of_range);
}

// -------------------------------------------------------------------------------------------- division

TEST(DivisionTest, SymmetricDivisionSplitsTheCellInTwoEqualHalvesWithoutMovingAnything)
{
	environment env = make_environment(3);
	const std::array<real_t, 3> centre = { 1.0, 2.0, 3.0 };
	auto* parent = env.create_cell(centre, 200);
	parent->is_motile() = 1;
	parent->migration_speed() = 2.5;
	parent->motility_direction()[1] = 1.0;

	const auto positions_before = positions_of(env);
	const std::vector<real_t> normal = { 1.0, 0.0, 0.0 };

	auto* daughter = env.divide_cell(0, normal);

	ASSERT_EQ(env.cells->size(), 2U);
	EXPECT_EQ(env.cells->get_agent_at(1), daughter);
	EXPECT_EQ(env.agents->size(), 200U);
	EXPECT_EQ(positions_of(env), positions_before);

	EXPECT_EQ(parent->agent_count(), index_t { 100 });
	EXPECT_EQ(daughter->agent_count(), index_t { 100 });

	// every sub-agent of the daughter lies at least as far along the normal as every one of the parent
	real_t highest_of_parent = std::numeric_limits<real_t>::lowest();
	real_t lowest_of_daughter = std::numeric_limits<real_t>::max();
	for (index_t i = 0; i < env.agents->size(); ++i)
	{
		const real_t p = projection(env, i, normal);
		if (env.agents->get_agent_at(i)->cell_index() == 0)
			highest_of_parent = std::max(highest_of_parent, p);
		else
			lowest_of_daughter = std::min(lowest_of_daughter, p);
	}
	EXPECT_LE(highest_of_parent, lowest_of_daughter);

	// the cells sit at the centres of their halves, on either side of the centre of the whole cell
	real_t centre_of_all = 0.0;
	for (index_t i = 0; i < 200; ++i)
	{
		centre_of_all += positions_before[3 * i] / 200.0;
	}
	EXPECT_LT(parent->position()[0], centre_of_all);
	EXPECT_GT(daughter->position()[0], centre_of_all);
	EXPECT_NEAR(0.5 * (parent->position()[0] + daughter->position()[0]), centre_of_all, 1e-9);
	EXPECT_NEAR(parent->position()[1], 2.0, 0.3);
	EXPECT_NEAR(daughter->position()[2], 3.0, 0.3);

	// the daughter inherits the properties of the parent
	EXPECT_EQ(daughter->is_motile(), 1);
	EXPECT_DOUBLE_EQ(daughter->migration_speed(), 2.5);
	EXPECT_DOUBLE_EQ(daughter->motility_direction()[1], 1.0);
}

class AsymmetricDivisionTest : public ::testing::TestWithParam<real_t>
{};

INSTANTIATE_TEST_SUITE_P(DaughterFractions, AsymmetricDivisionTest, ::testing::Values(0.1, 0.3, 0.5, 0.7, 0.9));

TEST_P(AsymmetricDivisionTest, DaughterTakesTheRequestedFractionFromTheSideTheNormalPointsTo)
{
	const real_t fraction = GetParam();
	environment env = make_environment(3);
	const std::array<real_t, 3> centre = { 0.0, 0.0, 0.0 };
	auto* parent = env.create_cell(centre, 200);

	// an oblique normal that is not normalised
	const std::vector<real_t> normal = { 1.0, 2.0, -0.5 };
	auto* daughter = env.divide_cell(0, normal, fraction);

	const auto expected = static_cast<index_t>(std::llround(fraction * 200.0));
	EXPECT_EQ(daughter->agent_count(), expected);
	EXPECT_EQ(parent->agent_count(), 200 - expected);

	// the daughter is the part further along the normal
	real_t along_normal = 0.0;
	for (std::size_t d = 0; d < 3; ++d)
	{
		along_normal += (daughter->position()[d] - parent->position()[d]) * normal[d];
	}
	EXPECT_GT(along_normal, 0.0);

	real_t highest_of_parent = std::numeric_limits<real_t>::lowest();
	real_t lowest_of_daughter = std::numeric_limits<real_t>::max();
	for (index_t i = 0; i < env.agents->size(); ++i)
	{
		const real_t p = projection(env, i, normal);
		if (env.agents->get_agent_at(i)->cell_index() == 0)
			highest_of_parent = std::max(highest_of_parent, p);
		else
			lowest_of_daughter = std::min(lowest_of_daughter, p);
	}
	EXPECT_LE(highest_of_parent, lowest_of_daughter);
}

TEST(DivisionTest, ReversingTheNormalSwapsWhichPartBecomesTheDaughter)
{
	const auto daughter_x = [](real_t normal_x) {
		environment env = make_environment(3);
		const std::array<real_t, 3> centre = { 0.0, 0.0, 0.0 };
		// 201 sub-agents are complete shells of the lattice, so the cell is the same seen from both sides
		env.create_cell(centre, 201);
		const std::vector<real_t> normal = { normal_x, 0.0, 0.0 };
		return env.divide_cell(0, normal, 0.3)->position()[0];
	};

	EXPECT_GT(daughter_x(1.0), 1.0);
	EXPECT_LT(daughter_x(-1.0), -1.0);
	// sub-agents exactly on the plane are shared out in container order, hence not exactly mirrored
	EXPECT_NEAR(daughter_x(1.0), -daughter_x(-1.0), 0.2);
}

TEST(DivisionTest, EachDaughterKeepsAtLeastOneSubAgent)
{
	environment env = make_environment(3);
	const std::array<real_t, 3> centre = { 0.0, 0.0, 0.0 };
	auto* parent = env.create_cell(centre, 5);
	const std::vector<real_t> normal = { 0.0, 0.0, 1.0 };

	// one percent of five sub-agents rounds to none
	auto* daughter = env.divide_cell(0, normal, 0.01);

	EXPECT_EQ(daughter->agent_count(), index_t { 1 });
	EXPECT_EQ(parent->agent_count(), index_t { 4 });
}

TEST(DivisionTest, OtherCellsAreNotTouched)
{
	environment env = make_environment(3);
	const std::array<real_t, 3> first = { 0.0, 0.0, 0.0 };
	const std::array<real_t, 3> second = { 20.0, 0.0, 0.0 };
	env.create_cell(first, 30);
	env.create_cell(second, 40);
	const std::vector<real_t> normal = { 0.0, 1.0, 0.0 };

	env.divide_cell(1, normal, 0.25);

	ASSERT_EQ(env.cells->size(), 3U);
	EXPECT_EQ(env.cells->get_agent_at(0)->agent_count(), index_t { 30 });
	EXPECT_EQ(env.cells->get_agent_at(1)->agent_count(), index_t { 30 });
	EXPECT_EQ(env.cells->get_agent_at(2)->agent_count(), index_t { 10 });

	const auto indices = cell_indices_of(env);
	EXPECT_EQ(std::ranges::count(indices, index_t { 0 }), 30);
	EXPECT_EQ(std::ranges::count(indices, index_t { 1 }), 30);
	EXPECT_EQ(std::ranges::count(indices, index_t { 2 }), 10);
	// the first thirty sub-agents belong to the first cell and stay there
	EXPECT_TRUE(std::ranges::all_of(indices.begin(), indices.begin() + 30, [](index_t c) { return c == 0; }));
}

TEST(DivisionTest, WorksInTwoDimensions)
{
	environment env = make_environment(2);
	const std::array<real_t, 2> centre = { 0.0, 0.0 };
	auto* parent = env.create_cell(centre, 60);
	const std::vector<real_t> normal = { 0.0, 1.0 };

	auto* daughter = env.divide_cell(0, normal, 0.3);

	EXPECT_EQ(daughter->agent_count(), index_t { 18 });
	EXPECT_EQ(parent->agent_count(), index_t { 42 });
	EXPECT_GT(daughter->position()[1], parent->position()[1]);
}

TEST(DivisionTest, RejectsInvalidArguments)
{
	environment env = make_environment(3);
	const std::array<real_t, 3> centre = { 0.0, 0.0, 0.0 };
	env.create_cell(centre, 10);
	env.create_cell(centre, 1);

	const std::vector<real_t> normal = { 1.0, 0.0, 0.0 };
	const std::vector<real_t> flat_normal = { 1.0, 0.0 };
	const std::vector<real_t> zero_normal = { 0.0, 0.0, 0.0 };

	EXPECT_THROW(env.divide_cell(5, normal), std::out_of_range);
	EXPECT_THROW(env.divide_cell(0, flat_normal), std::invalid_argument);
	EXPECT_THROW(env.divide_cell(0, zero_normal), std::invalid_argument);
	EXPECT_THROW(env.divide_cell(0, normal, 0.0), std::invalid_argument);
	EXPECT_THROW(env.divide_cell(0, normal, 1.0), std::invalid_argument);
	EXPECT_THROW(env.divide_cell(0, normal, -0.2), std::invalid_argument);
	EXPECT_THROW(env.divide_cell(0, normal, std::numeric_limits<real_t>::quiet_NaN()), std::invalid_argument);
	// a single sub-agent cannot be shared between two cells
	EXPECT_THROW(env.divide_cell(1, normal), std::invalid_argument);

	// nothing was changed by the failed calls
	EXPECT_EQ(env.cells->size(), 2U);
	EXPECT_EQ(env.cells->get_agent_at(0)->agent_count(), index_t { 10 });
}

// ---------------------------------------------------------------------------------------- longest axis

TEST(LongestAxisTest, FindsTheDirectionAStretchedCellIsElongatedAlong)
{
	environment env = make_environment(3);
	env.create_cell();

	// a rod along the direction (3, 4, 0) / 5 with a little thickness
	for (int n = -10; n <= 10; ++n)
	{
		for (const real_t side : { -0.2, 0.2 })
		{
			auto position = env.create_agent(0)->position();
			position[0] = 0.6 * n - 0.8 * side;
			position[1] = 0.8 * n + 0.6 * side;
			position[2] = side;
		}
	}

	const auto axis = env.longest_axis(0);

	ASSERT_EQ(axis.size(), 3U);
	EXPECT_NEAR(axis[0], 0.6, 1e-6);
	EXPECT_NEAR(axis[1], 0.8, 1e-6);
	EXPECT_NEAR(axis[2], 0.0, 1e-6);
}

TEST(LongestAxisTest, IsAUnitVectorWithItsLargestComponentPositive)
{
	environment env = make_environment(3);
	env.create_cell();
	// a rod along -z: an axis has no sign, the convention makes it +z
	for (int n = 0; n < 8; ++n)
	{
		env.create_agent(0)->position()[2] = -1.0 * n;
	}

	const auto axis = env.longest_axis(0);

	EXPECT_NEAR(axis[0], 0.0, 1e-12);
	EXPECT_NEAR(axis[1], 0.0, 1e-12);
	EXPECT_NEAR(axis[2], 1.0, 1e-12);
}

TEST(LongestAxisTest, DividingAcrossTheLongestAxisGivesTwoCompactDaughters)
{
	environment env = make_environment(2);
	auto* parent = env.create_cell();
	// a 12 x 4 rectangle of sub-agents
	for (int x = 0; x < 12; ++x)
	{
		for (int y = 0; y < 4; ++y)
		{
			auto position = env.create_agent(0)->position();
			position[0] = x;
			position[1] = y;
		}
	}

	const auto axis = env.longest_axis(0);
	EXPECT_NEAR(axis[0], 1.0, 1e-9);
	EXPECT_NEAR(axis[1], 0.0, 1e-9);

	auto* daughter = env.divide_cell(0, axis);

	// two 6 x 4 halves side by side
	EXPECT_EQ(parent->agent_count(), index_t { 24 });
	EXPECT_EQ(daughter->agent_count(), index_t { 24 });
	EXPECT_NEAR(parent->position()[0], 2.5, 1e-9);
	EXPECT_NEAR(daughter->position()[0], 8.5, 1e-9);
	EXPECT_NEAR(parent->position()[1], 1.5, 1e-9);
	EXPECT_NEAR(daughter->position()[1], 1.5, 1e-9);
}

TEST(LongestAxisTest, DegenerateCellsStillGiveAUnitVector)
{
	environment env = make_environment(3);
	env.create_cell();
	env.create_cell();
	env.create_agent(1);

	for (const index_t cell_index : { 0, 1 })
	{
		const auto axis = env.longest_axis(cell_index);
		real_t squared_length = 0.0;
		for (const real_t x : axis)
		{
			squared_length += x * x;
		}
		EXPECT_NEAR(squared_length, 1.0, 1e-12);
	}
	EXPECT_THROW(env.longest_axis(2), std::out_of_range);
}

} // namespace physicore::mechanics::micromechanics::tests
