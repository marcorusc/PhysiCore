#include <algorithm>
#include <array>

#include <gtest/gtest.h>

#include "cell_solver.h"
#include "openmp_solver.h"
#include "test_helpers.h"

namespace physicore::mechanics::micromechanics::kernels::openmp_solver::tests {

namespace {
void set_velocity(agent_interface* element, std::initializer_list<real_t> velocity)
{
	std::ranges::copy(velocity, element->velocity().begin());
}
} // namespace

TEST(CellSolverTest, PositionAndVelocityAreMeansOverTheSubAgents)
{
	environment env = make_environment(3);
	auto* first = env.create_cell();
	auto* second = env.create_cell();

	// created interleaved, so the cells do not own contiguous ranges of sub-agents
	set_velocity(add_agent(env, 0, { 1.0, 2.0, 3.0 }), { 1.0, 0.0, 0.0 });
	set_velocity(add_agent(env, 1, { -4.0, 0.0, 8.0 }), { 0.0, 6.0, 0.0 });
	set_velocity(add_agent(env, 0, { 3.0, 2.0, -1.0 }), { 3.0, 0.0, -2.0 });
	set_velocity(add_agent(env, 0, { 5.0, -1.0, 4.0 }), { 2.0, 3.0, -1.0 });

	cell_solver solver;
	solver.update(env);

	const std::array<real_t, 3> expected_position = { 3.0, 1.0, 2.0 };
	const std::array<real_t, 3> expected_velocity = { 2.0, 1.0, -1.0 };
	for (std::size_t d = 0; d < 3; ++d)
	{
		EXPECT_DOUBLE_EQ(first->position()[d], expected_position.at(d));
		EXPECT_DOUBLE_EQ(first->velocity()[d], expected_velocity.at(d));
	}
	EXPECT_EQ(first->agent_count(), index_t { 3 });

	// a cell with a single sub-agent coincides with it
	EXPECT_DOUBLE_EQ(second->position()[0], -4.0);
	EXPECT_DOUBLE_EQ(second->position()[2], 8.0);
	EXPECT_DOUBLE_EQ(second->velocity()[1], 6.0);
	EXPECT_EQ(second->agent_count(), index_t { 1 });
}

TEST(CellSolverTest, WorksInTwoDimensions)
{
	environment env = make_environment(2);
	auto* cell = env.create_cell();
	add_agent(env, 0, { 0.0, 4.0 });
	add_agent(env, 0, { 2.0, -2.0 });

	cell_solver solver;
	solver.update(env);

	EXPECT_DOUBLE_EQ(cell->position()[0], 1.0);
	EXPECT_DOUBLE_EQ(cell->position()[1], 1.0);
}

TEST(CellSolverTest, CellWithoutSubAgentsKeepsItsValues)
{
	environment env = make_environment(3);
	auto* empty = env.create_cell();
	env.create_cell();
	empty->position()[0] = 7.0;
	empty->velocity()[1] = -2.0;
	add_agent(env, 1, { 1.0, 1.0, 1.0 });

	cell_solver solver;
	solver.update(env);

	EXPECT_EQ(empty->agent_count(), index_t { 0 });
	EXPECT_DOUBLE_EQ(empty->position()[0], 7.0);
	EXPECT_DOUBLE_EQ(empty->velocity()[1], -2.0);
}

TEST(CellSolverTest, SubAgentsWithoutACellAreIgnored)
{
	environment env = make_environment(3);
	auto* cell = env.create_cell();
	add_agent(env, 0, { 2.0, 0.0, 0.0 });
	env.agents->create()->position()[0] = 100.0;

	cell_solver solver;
	solver.update(env);

	EXPECT_EQ(cell->agent_count(), index_t { 1 });
	EXPECT_DOUBLE_EQ(cell->position()[0], 2.0);
}

TEST(CellSolverTest, SolverKeepsCellsInSyncWithTheirSubAgents)
{
	environment env = make_environment(3, 0.002);
	add_jittered_lattice_cell(env, 1.0, 0.15, 0.0, 42);
	add_jittered_lattice_cell(env, 1.0, 0.15, 3.3, 43);
	auto* left = env.cells->get_agent_at(0);
	auto* right = env.cells->get_agent_at(1);

	openmp_solver solver;
	solver.initialize(env);

	// initialize alone already places the cells at the centres of their lattices
	EXPECT_NEAR(left->position()[0], 1.0, 0.1);
	EXPECT_NEAR(right->position()[0], 4.3, 0.1);
	const real_t initial_distance = right->position()[0] - left->position()[0];

	solver.solve(env, 500);

	// the cells attract each other: they move towards one another and their velocities say so
	EXPECT_LT(right->position()[0] - left->position()[0], initial_distance - 0.05);
	EXPECT_GT(left->velocity()[0], 0.0);
	EXPECT_LT(right->velocity()[0], 0.0);

	// and the centre of mass of the pair stays where it was, as there is no external force
	const auto centre = centre_of_mass(env);
	for (std::size_t d = 0; d < 3; ++d)
	{
		EXPECT_NEAR(0.5 * (left->position()[d] + right->position()[d]), centre.at(d), 1e-10);
	}
}

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver::tests
