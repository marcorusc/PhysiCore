#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "openmp_solver.h"
#include "test_helpers.h"

// Growth and division together with the solver: a cell that gains elements must stay one compact cell
// of the right size, and the two parts of a divided cell must stay two cells.

namespace physicore::mechanics::micromechanics::kernels::openmp_solver::tests {

namespace {
struct cell_shape
{
	index_t count = 0;
	std::array<real_t, 3> centre {};
	// sqrt(5/3) times the radius of gyration: the radius of a uniform ball with the same spread
	real_t radius = 0.0;
	// largest distance of a sub-agent from the centre
	real_t furthest = 0.0;
};

cell_shape shape_of(environment& env, index_t cell_index)
{
	std::vector<index_t> members;
	for (index_t i = 0; i < env.agents->size(); ++i)
	{
		if (env.agents->get_agent_at(i)->cell_index() == cell_index)
		{
			members.push_back(i);
		}
	}

	cell_shape shape;
	shape.count = members.size();
	for (const index_t i : members)
	{
		const auto position = env.agents->get_agent_at(i)->position();
		for (std::size_t d = 0; d < 3; ++d)
		{
			shape.centre.at(d) += position[d] / static_cast<real_t>(members.size());
		}
	}

	real_t squared_gyration_radius = 0.0;
	for (const index_t i : members)
	{
		const auto position = env.agents->get_agent_at(i)->position();
		real_t squared = 0.0;
		for (std::size_t d = 0; d < 3; ++d)
		{
			squared += std::pow(position[d] - shape.centre.at(d), 2);
		}
		squared_gyration_radius += squared / static_cast<real_t>(members.size());
		shape.furthest = std::max(shape.furthest, std::sqrt(squared));
	}
	shape.radius = std::sqrt(5.0 / 3.0 * squared_gyration_radius);

	return shape;
}

// Largest distance from any sub-agent to its nearest neighbour: small when nobody has been thrown out
real_t loneliest(environment& env)
{
	real_t result = 0.0;
	for (index_t i = 0; i < env.agents->size(); ++i)
	{
		real_t nearest = std::numeric_limits<real_t>::max();
		for (index_t j = 0; j < env.agents->size(); ++j)
		{
			if (i != j)
			{
				nearest = std::min(nearest, squared_distance(env, i, j));
			}
		}
		result = std::max(result, std::sqrt(nearest));
	}
	return result;
}

// Timestep at a twentieth of the stability limit of an isolated pair, for unit well depth and drag
real_t timestep_for(real_t steepness) { return 0.05 / (8.0 * steepness * steepness); }

index_t steps_for(real_t time, real_t steepness) { return static_cast<index_t>(time / timestep_for(steepness)); }

// Adds sub-agents to a cell one at a time, letting the cell relax briefly after each
void grow(environment& env, openmp_solver& solver, index_t cell_index, index_t additions, real_t steepness)
{
	for (index_t n = 0; n < additions; ++n)
	{
		env.grow_cell(cell_index);
		solver.solve(env, steps_for(0.05, steepness));
	}
}

// Separation of two cells along a direction: smallest projection of the second minus largest of the first
real_t gap_along(environment& env, const std::vector<real_t>& direction, index_t first, index_t second)
{
	real_t highest_of_first = std::numeric_limits<real_t>::lowest();
	real_t lowest_of_second = std::numeric_limits<real_t>::max();
	for (index_t i = 0; i < env.agents->size(); ++i)
	{
		auto* element = env.agents->get_agent_at(i);
		real_t projection = 0.0;
		for (std::size_t d = 0; d < 3; ++d)
		{
			projection += element->position()[d] * direction[d];
		}
		if (element->cell_index() == first)
			highest_of_first = std::max(highest_of_first, projection);
		else if (element->cell_index() == second)
			lowest_of_second = std::min(lowest_of_second, projection);
	}
	return lowest_of_second - highest_of_first;
}

real_t distance_between(const cell_shape& a, const cell_shape& b)
{
	real_t squared = 0.0;
	for (std::size_t d = 0; d < 3; ++d)
	{
		squared += std::pow(a.centre.at(d) - b.centre.at(d), 2);
	}
	return std::sqrt(squared);
}
} // namespace

class GrowthWithSolverTest : public ::testing::TestWithParam<real_t>
{};

INSTANTIATE_TEST_SUITE_P(Steepness, GrowthWithSolverTest, ::testing::Values(2.0, 4.0));

TEST_P(GrowthWithSolverTest, DoublingTheSubAgentsGrowsTheCellWithoutThrowingAnyOut)
{
	const real_t steepness = GetParam();

	environment env = make_environment(3, timestep_for(steepness));
	env.params.steepness = steepness;
	const std::array<real_t, 3> centre = { 0.0, 0.0, 0.0 };
	env.create_cell(centre, 32);

	openmp_solver solver;
	solver.solve(env, steps_for(0.5, steepness));
	const cell_shape before = shape_of(env, 0);

	grow(env, solver, 0, 32, steepness);
	solver.solve(env, steps_for(1.0, steepness));
	const cell_shape after = shape_of(env, 0);


	EXPECT_EQ(after.count, index_t { 64 });
	EXPECT_EQ(env.cells->get_agent_at(0)->agent_count(), index_t { 64 });

	// every sub-agent still has a neighbour about one equilibrium distance away
	EXPECT_LT(loneliest(env), 1.2);

	// Twice the sub-agents is twice the volume, so a compact cell grows in radius by 2^(1/3) = 1.26.
	// New sub-agents attach where there is room, which is at the surface, and without noise only a
	// soft potential lets them settle inwards. Measured: 1.35 for steepness 2 with this short
	// relaxation (1.27 with a longer one) and 1.60 for steepness 4, whose cell stays lumpy.
	const real_t ratio = after.radius / before.radius;
	EXPECT_GT(ratio, 1.15);
	EXPECT_LT(ratio, steepness < 3.0 ? 1.5 : 2.0);
}

class DivisionWithSolverTest : public ::testing::TestWithParam<real_t>
{};

INSTANTIATE_TEST_SUITE_P(DaughterFractions, DivisionWithSolverTest, ::testing::Values(0.5, 0.3));

TEST_P(DivisionWithSolverTest, DaughtersStayTwoCompactCellsThatDoNotMix)
{
	const real_t fraction = GetParam();
	const real_t steepness = 2.0;

	environment env = make_environment(3, timestep_for(steepness));
	const std::array<real_t, 3> centre = { 0.0, 0.0, 0.0 };
	env.create_cell(centre, 80);

	openmp_solver solver;
	solver.solve(env, steps_for(1.0, steepness));

	const auto axis = env.longest_axis(0);
	env.divide_cell(0, axis, fraction);

	const real_t distance_at_division = distance_between(shape_of(env, 0), shape_of(env, 1));
	const real_t gap_at_division = gap_along(env, axis, 0, 1);

	solver.solve(env, steps_for(6.0, steepness));

	const cell_shape parent = shape_of(env, 0);
	const cell_shape daughter = shape_of(env, 1);
	const real_t gap = gap_along(env, axis, 0, 1);


	const auto expected_daughter = static_cast<index_t>(std::llround(fraction * 80.0));
	EXPECT_EQ(daughter.count, expected_daughter);
	EXPECT_EQ(parent.count, 80 - expected_daughter);
	EXPECT_EQ(env.cells->get_agent_at(0)->agent_count(), parent.count);
	EXPECT_EQ(env.cells->get_agent_at(1)->agent_count(), daughter.count);

	// the plane separates the daughters cleanly at the moment of division
	EXPECT_GE(gap_at_division, 0.0);
	// afterwards their surfaces may interlock, but by less than one element (measured: +0.35 for the
	// symmetric division, -0.36 for the asymmetric one)
	EXPECT_GT(gap, -1.0);
	// the daughters do not collapse back into each other
	EXPECT_GT(distance_between(parent, daughter), 0.9 * distance_at_division);
	// each daughter is a compact cell and nobody was thrown out
	EXPECT_LT(parent.furthest, 1.5 * parent.radius);
	EXPECT_LT(daughter.furthest, 1.5 * daughter.radius);
	EXPECT_LT(loneliest(env), 1.2);
}

TEST(GenerationsTest, RepeatedGrowthAndDivisionIsReproducible)
{
	// One cell of 16 sub-agents grows to 32 and divides, twice over: four cells of 16.
	const auto run = [](std::uint64_t seed) {
		const real_t steepness = 2.0;
		environment env = make_environment(3, timestep_for(steepness));
		env.params.random_seed = seed;
		const std::array<real_t, 3> centre = { 0.0, 0.0, 0.0 };
		env.create_cell(centre, 16);

		openmp_solver solver;
		for (int generation = 0; generation < 2; ++generation)
		{
			const index_t cells_count = env.cells->size();
			for (index_t c = 0; c < cells_count; ++c)
			{
				grow(env, solver, c, 16, steepness);
			}
			for (index_t c = 0; c < cells_count; ++c)
			{
				env.divide_cell(c, env.longest_axis(c));
			}
			solver.solve(env, steps_for(1.0, steepness));
		}

		std::vector<index_t> counts;
		counts.reserve(env.cells->size());
		for (index_t c = 0; c < env.cells->size(); ++c)
		{
			counts.push_back(env.cells->get_agent_at(c)->agent_count());
		}
		return std::pair { counts, positions_of(env) };
	};

	const auto [counts, positions] = run(11);

	EXPECT_EQ(counts, std::vector<index_t>({ 16, 16, 16, 16 }));
	EXPECT_EQ(positions.size(), 64U * 3U);

	const auto [same_counts, same_positions] = run(11);
	EXPECT_EQ(same_positions, positions);

	const auto [other_counts, other_positions] = run(12);
	EXPECT_EQ(other_counts, counts);
	EXPECT_NE(other_positions, positions);
}

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver::tests
