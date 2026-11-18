#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>
#include <micromechanics/cell_type.h>

#include "force_solver.h"
#include "openmp_solver.h"
#include "test_helpers.h"

// Cells of different types differ in the size, stiffness and drag of their elements and in how
// strongly they stick to each other. These tests check that the solver uses the right interaction
// for every pair of elements. What cell types mean for noise, motility, pressure and walls is
// tested with those features.

namespace physicore::mechanics::micromechanics::kernels::openmp_solver::tests {

// ------------------------------------------------------------------------- which potential acts

TEST(CellTypesTest, ElementsOfDifferentTypesRestAtTheSumOfTheirRadii)
{
	environment env = make_environment(3);
	const index_t large = env.add_cell_type(large_type());
	const interaction_table table(env.cell_types, env.params);

	// radii 0.5 and 1: neither attraction nor repulsion at the distance 1.5
	for (const real_t distance : { 1.2, 1.5, 1.8 })
	{
		environment trial = make_environment(3);
		trial.add_cell_type(large_type());
		add_point_cell(trial, 0, 0.0);
		add_point_cell(trial, large, distance);

		const auto forces = forces_of(trial);
		const real_t expected = force_coefficient(table.inter(0, large), distance * distance) * distance;

		// on the second cell, along +x: positive pushes it away
		EXPECT_DOUBLE_EQ(forces[3], expected) << "distance " << distance;
		// action equals reaction
		EXPECT_DOUBLE_EQ(forces[0], -forces[3]);
	}

	environment at_rest = make_environment(3);
	at_rest.add_cell_type(large_type());
	add_point_cell(at_rest, 0, 0.0);
	add_point_cell(at_rest, large, 1.5);
	EXPECT_NEAR(forces_of(at_rest)[3], 0.0, 1e-14);
}

TEST(CellTypesTest, ThreeKindsOfPairsUseThreeDifferentWells)
{
	// one cell of type 0 with two elements, a second cell of type 0 and a cell of another type, all
	// at the same distance of 1.2 from the first element
	environment env = make_environment(3);
	cell_type other;
	other.inter_well_depth = 2.0;
	other.adhesion_affinities = { 0.5, 1.0 };
	const index_t other_index = env.add_cell_type(other);
	env.cell_types[0].adhesion_affinities = { 0.8, 0.5 };

	const interaction_table table(env.cell_types, env.params);
	const real_t distance = 1.2;
	const real_t squared = distance * distance;

	const auto force_on_first = [&](auto&& add_partner) {
		environment trial = make_environment(3);
		trial.cell_types = env.cell_types;
		trial.create_cell();
		add_agent(trial, 0, { 0.0, 0.0, 0.0 });
		add_partner(trial);
		return forces_of(trial)[0];
	};

	// same cell: the intra-cellular well, in full
	const real_t same_cell = force_on_first([&](environment& trial) { add_agent(trial, 0, { distance, 0.0, 0.0 }); });
	EXPECT_DOUBLE_EQ(same_cell, -force_coefficient(intra_potential(env.cell_types[0], env.params), squared) * distance);

	// another cell of the same type: the inter-cellular well, attraction times the affinity 0.8
	const real_t same_type = force_on_first([&](environment& trial) { add_point_cell(trial, 0, distance); });
	EXPECT_DOUBLE_EQ(same_type,
					 -0.8 * force_coefficient(inter_potential(env.cell_types[0], env.params), squared) * distance);

	// a cell of the other type: the mixed well
	const real_t other_type = force_on_first([&](environment& trial) { add_point_cell(trial, other_index, distance); });
	EXPECT_DOUBLE_EQ(other_type, -force_coefficient(table.inter(0, other_index), squared) * distance);

	// all three attract (the first element is pulled towards +x) with different strengths
	EXPECT_GT(same_cell, 0.0);
	EXPECT_GT(same_type, 0.0);
	EXPECT_GT(other_type, 0.0);
	EXPECT_NE(same_cell, same_type);
	EXPECT_NE(same_type, other_type);
}

TEST(CellTypesTest, ForcesBetweenCellsOfDifferentTypesAddUpToZero)
{
	environment env = make_environment(3);
	cell_type large = large_type();
	large.adhesion_affinities = { 0.3, 1.0 };
	const index_t large_index = env.add_cell_type(large);
	env.cell_types[0].adhesion_affinities = { 1.0, 0.7 };

	const std::array<real_t, 3> first_centre = { 0.0, 0.0, 0.0 };
	const std::array<real_t, 3> second_centre = { 3.2, 0.4, -0.3 };
	env.create_cell_of_type(0, first_centre, 13);
	env.create_cell_of_type(large_index, second_centre, 13);

	const auto forces = forces_of(env);

	std::array<real_t, 3> total {};
	real_t largest = 0.0;
	for (std::size_t k = 0; k < forces.size(); ++k)
	{
		total.at(k % 3) += forces[k];
		largest = std::max(largest, std::abs(forces[k]));
	}
	ASSERT_GT(largest, 0.1);
	for (const real_t component : total)
	{
		EXPECT_NEAR(component, 0.0, 1e-12 * largest * static_cast<real_t>(forces.size()));
	}
}

TEST(CellTypesTest, SubAgentWithoutACellInteractsLikeTypeZero)
{
	const auto force_on_loose = [](bool as_cell) {
		environment env = make_environment(3);
		const index_t large = env.add_cell_type(large_type());
		add_point_cell(env, large, 0.0);
		if (as_cell)
			add_point_cell(env, 0, 1.3);
		else
			env.agents->create()->position()[0] = 1.3;
		return forces_of(env)[3];
	};

	EXPECT_DOUBLE_EQ(force_on_loose(false), force_on_loose(true));
	EXPECT_NE(force_on_loose(false), 0.0);
}

// ------------------------------------------------------------------------------------- affinity

TEST(CellTypesTest, AffinityZeroLeavesRepulsionButNoAttraction)
{
	const auto final_distance = [](real_t affinity, real_t start) {
		environment env = make_environment(3, 0.002);
		env.cell_types[0].adhesion_affinities = { affinity };
		add_point_cell(env, 0, 0.0);
		add_point_cell(env, 0, start);

		openmp_solver solver;
		solver.initialize(env);
		solver.solve(env, 20000);
		return x_of(env, 1) - x_of(env, 0);
	};

	// started apart: full affinity pulls the two together, zero affinity leaves them where they are
	EXPECT_NEAR(final_distance(1.0, 1.4), 1.0, 1e-6);
	EXPECT_DOUBLE_EQ(final_distance(0.0, 1.4), 1.4);

	// started overlapping: both are pushed apart to the equilibrium distance, and no further
	EXPECT_NEAR(final_distance(1.0, 0.6), 1.0, 1e-6);
	EXPECT_NEAR(final_distance(0.0, 0.6), 1.0, 1e-3);
	EXPECT_LE(final_distance(0.0, 0.6), 1.0 + 1e-12);
}

TEST(CellTypesTest, WeakerAffinityPullsCellsTogetherMoreSlowly)
{
	const auto approach_after = [](real_t affinity) {
		environment env = make_environment(3, 0.002);
		env.cell_types[0].adhesion_affinities = { affinity };
		add_point_cell(env, 0, 0.0);
		add_point_cell(env, 0, 1.4);

		openmp_solver solver;
		solver.initialize(env);
		solver.solve(env, 100);
		return 1.4 - (x_of(env, 1) - x_of(env, 0));
	};

	const real_t full = approach_after(1.0);
	const real_t half = approach_after(0.5);

	EXPECT_GT(half, 0.0);
	EXPECT_GT(full, half);
	// in the first step the speed is proportional to the affinity
	environment env = make_environment(3, 0.002);
	env.cell_types[0].adhesion_affinities = { 0.5 };
	add_point_cell(env, 0, 0.0);
	add_point_cell(env, 0, 1.4);
	const real_t half_force = forces_of(env)[3];
	env.cell_types[0].adhesion_affinities = { 1.0 };
	EXPECT_DOUBLE_EQ(half_force, 0.5 * forces_of(env)[3]);
}

// ------------------------------------------------------------------------------------------ drag

TEST(CellTypesTest, VelocityOfAnElementIsItsForceOverTheDragOfItsType)
{
	environment env = make_environment(3, 0.001);
	cell_type heavy;
	heavy.drag_coefficient = 4.0;
	const index_t heavy_index = env.add_cell_type(heavy);

	add_point_cell(env, 0, 0.0);
	add_point_cell(env, heavy_index, 0.8); // overlapping: they repel with equal and opposite forces

	const real_t force = forces_of(env)[3];
	ASSERT_GT(force, 0.0);

	openmp_solver solver;
	solver.initialize(env);
	solver.solve(env, 1);

	EXPECT_DOUBLE_EQ(env.agents->get_agent_at(0)->velocity()[0], -force / 1.0);
	EXPECT_DOUBLE_EQ(env.agents->get_agent_at(1)->velocity()[0], force / 4.0);
}

// -------------------------------------------------------------------------------------- checking

TEST(CellTypesTest, SolverRefusesACellOfATypeThatDoesNotExistAndTypesThatMakeNoSense)
{
	environment env = make_environment(3);
	env.create_cell();
	openmp_solver solver;
	EXPECT_NO_THROW(solver.initialize(env));

	env.cells->get_agent_at(0)->type_index() = 3;
	EXPECT_THROW(solver.initialize(env), std::invalid_argument);

	env.cells->get_agent_at(0)->type_index() = 0;
	env.cell_types[0].drag_coefficient = 0.0;
	EXPECT_THROW(solver.initialize(env), std::invalid_argument);
}

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver::tests
