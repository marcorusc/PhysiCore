#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>
#include <micromechanics/solver_registry.h>

#include "openmp_solver.h"
#include "test_helpers.h"

namespace physicore::mechanics::micromechanics::kernels::openmp_solver::tests {

TEST(OpenMPSolverTest, SingleStepIsForwardEuler)
{
	environment env = make_environment(3, 0.01);
	env.cell_types[0].drag_coefficient = 2.0;
	env.create_cell();
	auto* left = add_agent(env, 0, { 0.0, 0.0, 0.0 });
	auto* right = add_agent(env, 0, { 0.8, 0.0, 0.0 });

	// repulsive: the pair is closer than the equilibrium distance of 1
	const real_t force_on_right = 0.8 * force_coefficient(intra_potential(env.cell_types[0], env.params), 0.8 * 0.8);
	ASSERT_GT(force_on_right, 0.0);

	openmp_solver solver;
	solver.initialize(env);
	solver.solve(env, 1);

	const real_t expected_velocity = force_on_right / 2.0;
	EXPECT_DOUBLE_EQ(right->velocity()[0], expected_velocity);
	EXPECT_DOUBLE_EQ(left->velocity()[0], -expected_velocity);
	EXPECT_DOUBLE_EQ(right->position()[0], 0.8 + 0.01 * expected_velocity);
	EXPECT_DOUBLE_EQ(left->position()[0], -0.01 * expected_velocity);
	for (const std::size_t d : { 1U, 2U })
	{
		EXPECT_DOUBLE_EQ(right->position()[d], 0.0);
		EXPECT_DOUBLE_EQ(right->velocity()[d], 0.0);
	}

	solver.solve(env, 1);

	// the pair has moved apart, so the repulsion is weaker
	EXPECT_LT(right->velocity()[0], expected_velocity);
	EXPECT_GT(right->velocity()[0], 0.0);
}

class PairRelaxationTest : public ::testing::TestWithParam<std::tuple<index_t, real_t>>
{};

INSTANTIATE_TEST_SUITE_P(DimensionsAndStartDistances, PairRelaxationTest,
						 ::testing::Combine(::testing::Values(2, 3), ::testing::Values(0.7, 1.3)));

TEST_P(PairRelaxationTest, PairRelaxesToTheEquilibriumDistance)
{
	const auto [dims, start_distance] = GetParam();

	// along the diagonal, so every axis takes part
	const real_t offset = start_distance / std::sqrt(static_cast<real_t>(dims));

	// stiffness k = 32 with the default parameters, so the timestep is far below eta / k = 0.03125
	environment env = make_environment(dims, 0.005);
	env.create_cell();
	if (dims == 2)
	{
		add_agent(env, 0, { 1.0, 2.0 });
		add_agent(env, 0, { 1.0 + offset, 2.0 + offset });
	}
	else
	{
		add_agent(env, 0, { 1.0, 2.0, 3.0 });
		add_agent(env, 0, { 1.0 + offset, 2.0 + offset, 3.0 + offset });
	}
	const auto initial_centre = centre_of_mass(env);

	openmp_solver solver;
	solver.solve(env, 4000);

	EXPECT_NEAR(std::sqrt(squared_distance(env, 0, 1)), env.cell_types[0].equilibrium_distance, 1e-10);

	const auto final_centre = centre_of_mass(env);
	for (std::size_t d = 0; d < 3; ++d)
	{
		EXPECT_NEAR(final_centre.at(d), initial_centre.at(d), 1e-12);
	}
	for (const real_t velocity : env.agents->get_agent_at(0)->velocity())
	{
		EXPECT_NEAR(velocity, 0.0, 1e-9);
	}
}

TEST(OpenMPSolverTest, EnergyNeverIncreases)
{
	// Overdamped dynamics without noise is a descent on the potential energy.
	environment env = make_environment(3, 0.001);
	add_jittered_lattice_cell(env, 1.0, 0.15);

	openmp_solver solver;

	const real_t initial_energy = total_energy(env);
	real_t previous_energy = initial_energy;
	for (index_t step = 0; step < 300; ++step)
	{
		solver.solve(env, 1);

		const real_t energy = total_energy(env);
		ASSERT_LE(energy, previous_energy + 1e-10) << "step " << step;
		previous_energy = energy;
	}

	EXPECT_LT(previous_energy, initial_energy - 0.1);
}

TEST(OpenMPSolverTest, CellStaysCohesiveAndKeepsItsCentreOfMass)
{
	environment env = make_environment(3, 0.002);
	add_jittered_lattice_cell(env, 1.0, 0.15);
	const index_t agents_count = env.agents->size();
	const auto initial_centre = centre_of_mass(env);

	openmp_solver solver;
	solver.solve(env, 3000);

	// no net force, so the centre of mass does not move
	const auto centre = centre_of_mass(env);
	for (std::size_t d = 0; d < 3; ++d)
	{
		EXPECT_NEAR(centre.at(d), initial_centre.at(d), 1e-10);
	}

	const auto positions = positions_of(env);
	for (index_t i = 0; i < agents_count; ++i)
	{
		// every element keeps a neighbour close to the equilibrium distance and overlaps with none
		real_t nearest = std::numeric_limits<real_t>::max();
		for (index_t j = 0; j < agents_count; ++j)
		{
			if (i != j)
			{
				nearest = std::min(nearest, std::sqrt(squared_distance(env, i, j)));
			}
		}
		EXPECT_GT(nearest, 0.8) << "element " << i;
		EXPECT_LT(nearest, 1.1) << "element " << i;

		// 27 close-packed elements of diameter 1 fill a ball of radius about 1.7
		real_t squared_radius = 0.0;
		for (std::size_t d = 0; d < 3; ++d)
		{
			squared_radius += std::pow(positions[3 * i + d] - centre.at(d), 2);
		}
		EXPECT_LT(std::sqrt(squared_radius), 2.2) << "element " << i;
	}

	// the cell has come to rest
	for (index_t i = 0; i < agents_count; ++i)
	{
		for (const real_t velocity : env.agents->get_agent_at(i)->velocity())
		{
			EXPECT_NEAR(velocity, 0.0, 1e-3);
		}
	}
}

TEST(OpenMPSolverTest, ImmovableCellKeepsItsPositionButStillPushesOthers)
{
	environment env = make_environment(3, 0.005);
	env.create_cell()->is_movable() = 0;
	env.create_cell();
	auto* fixed = add_agent(env, 0, { 0.0, 0.0, 0.0 });
	auto* free_agent = add_agent(env, 1, { 0.7, 0.0, 0.0 });

	openmp_solver solver;
	solver.solve(env, 4000);

	EXPECT_DOUBLE_EQ(fixed->position()[0], 0.0);
	EXPECT_DOUBLE_EQ(fixed->velocity()[0], 0.0);
	// only the free element moves, so it travels the whole way to the equilibrium distance
	EXPECT_NEAR(free_agent->position()[0], 1.0, 1e-9);
}

TEST(OpenMPSolverTest, SolvingSeveralIterationsEqualsRepeatedSingleSteps)
{
	environment at_once = make_environment(3, 0.002);
	environment step_by_step = make_environment(3, 0.002);
	add_jittered_lattice_cell(at_once, 1.0, 0.15);
	add_jittered_lattice_cell(step_by_step, 1.0, 0.15);

	openmp_solver first;
	openmp_solver second;
	first.solve(at_once, 25);
	for (index_t i = 0; i < 25; ++i)
	{
		second.solve(step_by_step, 1);
	}

	EXPECT_EQ(positions_of(at_once), positions_of(step_by_step));
}

TEST(OpenMPSolverTest, IsRegisteredAndRunsThroughTheEnvironment)
{
	environment env = make_environment(3, 0.005);
	env.create_cell();
	add_agent(env, 0, { 0.0, 0.0, 0.0 });
	add_agent(env, 0, { 0.7, 0.0, 0.0 });

	env.solver = solver_registry::instance().get("openmp_solver");
	ASSERT_NE(env.solver, nullptr);
	env.solver->initialize(env);

	env.run_single_timestep();

	EXPECT_GT(std::sqrt(squared_distance(env, 0, 1)), 0.7);
}

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver::tests
