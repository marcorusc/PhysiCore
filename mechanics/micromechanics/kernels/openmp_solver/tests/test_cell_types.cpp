#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>
#include <micromechanics/cell_type.h>
#include <micromechanics/domain_wall.h>

#include "cell_solver.h"
#include "force_solver.h"
#include "openmp_solver.h"
#include "pressure_solver.h"
#include "test_helpers.h"

// Cells of different types differ in the size, stiffness and drag of their elements and in how
// strongly they stick to each other. These tests check that the solver uses the right interaction
// for every pair of elements, and what follows from it for cells.

namespace physicore::mechanics::micromechanics::kernels::openmp_solver::tests {

namespace {
// Net force on every sub-agent, agents_count * dims values
std::vector<real_t> forces_of(environment& env)
{
	force_solver forces;
	forces.compute(env);
	return { forces.forces().begin(), forces.forces().end() };
}

// A cell of one sub-agent of the given type on the x axis
cell_interface* add_point_cell(environment& env, index_t type_index, real_t x)
{
	auto* cell = env.create_cell(type_index);
	env.create_agent(env.cells->size() - 1)->position()[0] = x;
	return cell;
}

real_t x_of(environment& env, index_t agent) { return env.agents->get_agent_at(agent)->position()[0]; }

// A second type whose elements are twice as large as those of the default type
cell_type large_type()
{
	cell_type type;
	type.name = "large";
	type.equilibrium_distance = 2.0;
	return type;
}
} // namespace

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

// ---------------------------------------------------------------------- one type is the old model

TEST(CellTypesTest, TypesWithIdenticalValuesBehaveLikeOneType)
{
	const auto run = [](bool two_types) {
		environment env = make_environment(3, 0.002);
		env.params.thermal_energy = 0.01;
		env.params.random_seed = 9;
		if (two_types)
		{
			// a second type with the same values. (A type with larger elements, even an unused one, would
			// widen the bins of the neighbour search and change the order of the sums, and with it the
			// last digits.)
			env.add_cell_type(cell_type {});
		}
		add_jittered_lattice_cell(env, 1.0, 0.1, 0.0, 42);
		add_jittered_lattice_cell(env, 1.0, 0.1, 3.3, 43);
		if (two_types)
		{
			env.cells->get_agent_at(1)->type_index() = 1;
		}

		openmp_solver solver;
		solver.initialize(env);
		solver.solve(env, 300);

		std::vector<real_t> result = positions_of(env);
		result.push_back(env.cells->get_agent_at(0)->pressure());
		return result;
	};

	// bitwise identical, not merely close
	EXPECT_EQ(run(true), run(false));
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

TEST(CellTypesTest, CellsWithoutAffinityStillCannotPassThroughEachOther)
{
	// Two cells on a line crawl head-on. With the affinity 0 they do not attract at all, but they
	// must still stop each other.
	const real_t speed = 0.2;
	environment env = make_environment(1, 0.002);
	env.cell_types[0].adhesion_affinities = { 0.0 };
	for (const real_t side : { -1.0, 1.0 })
	{
		const std::array<real_t, 1> centre = { 5.0 * side };
		auto* cell = env.create_cell(centre, 5);
		cell->is_motile() = 1;
		cell->migration_speed() = speed;
		cell->persistence_time() = 1.0;
		cell->migration_bias() = 1.0;
		cell->migration_bias_direction()[0] = -side;
	}

	openmp_solver solver;
	solver.initialize(env);
	solver.solve(env, 20000);

	real_t rightmost_of_left = -1e30;
	real_t leftmost_of_right = 1e30;
	for (index_t i = 0; i < env.agents->size(); ++i)
	{
		auto* element = env.agents->get_agent_at(i);
		if (element->cell_index() == 0)
			rightmost_of_left = std::max(rightmost_of_left, element->position()[0]);
		else
			leftmost_of_right = std::min(leftmost_of_right, element->position()[0]);
	}
	EXPECT_NEAR(env.cells->get_agent_at(0)->velocity()[0], 0.0, 0.01 * speed);
	EXPECT_GT(leftmost_of_right - rightmost_of_left, 0.5);
	EXPECT_LT(leftmost_of_right - rightmost_of_left, 1.0);

	// when they stop crawling they are pushed apart to the equilibrium distance and nothing holds
	// them together
	env.cells->get_agent_at(0)->is_motile() = 0;
	env.cells->get_agent_at(1)->is_motile() = 0;
	solver.solve(env, 20000);

	EXPECT_NEAR(env.cells->get_agent_at(0)->pressure(), 0.0, 1e-3);
}

// ---------------------------------------------------------------------------- stiffness and drag

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

TEST(CellTypesTest, SofterCellIsCompressedMoreThanAStifferOne)
{
	// A soft and a stiff cell on a line crawl against each other with the same traction. Both end up
	// shorter than at rest; the soft one by more.
	const real_t speed = 0.5;
	const index_t elements = 6;

	const auto lengths = [&](bool crawling) {
		environment env = make_environment(1, 0.001);
		env.cell_types[0].intra_well_depth = 0.5; // soft
		cell_type stiff;
		stiff.intra_well_depth = 4.0;
		const index_t stiff_index = env.add_cell_type(stiff);

		const std::array<index_t, 2> types = { 0, stiff_index };
		for (std::size_t k = 0; k < 2; ++k)
		{
			const real_t side = k == 0 ? -1.0 : 1.0;
			const std::array<real_t, 1> centre = { 6.0 * side };
			auto* cell = env.create_cell_of_type(types.at(k), centre, elements);
			cell->is_motile() = crawling ? 1 : 0;
			cell->migration_speed() = speed;
			cell->persistence_time() = 1.0;
			cell->migration_bias() = 1.0;
			cell->migration_bias_direction()[0] = -side;
		}

		openmp_solver solver;
		solver.initialize(env);
		solver.solve(env, 40000);

		std::array<real_t, 2> low = { 1e30, 1e30 };
		std::array<real_t, 2> high = { -1e30, -1e30 };
		for (index_t i = 0; i < env.agents->size(); ++i)
		{
			auto* element = env.agents->get_agent_at(i);
			low.at(element->cell_index()) = std::min(low.at(element->cell_index()), element->position()[0]);
			high.at(element->cell_index()) = std::max(high.at(element->cell_index()), element->position()[0]);
		}
		return std::array<real_t, 2> { high[0] - low[0], high[1] - low[1] };
	};

	const auto at_rest = lengths(false);
	const auto pressed = lengths(true);

	const real_t soft_compression = at_rest[0] - pressed[0];
	const real_t stiff_compression = at_rest[1] - pressed[1];

	EXPECT_GT(stiff_compression, 0.0);
	EXPECT_GT(soft_compression, 2.0 * stiff_compression);
}

TEST(CellTypesTest, EveryTypeDiffusesWithTheThermalEnergyOverItsDrag)
{
	const index_t per_type = 3000;
	const real_t thermal_energy = 0.2;
	const real_t timestep = 0.01;
	const index_t steps = 200;

	environment env = make_environment(3, timestep);
	env.params.thermal_energy = thermal_energy;
	env.params.random_seed = 21;
	env.cell_types[0].inter_well_depth = 0.0; // free elements
	cell_type heavy;
	heavy.drag_coefficient = 4.0;
	heavy.inter_well_depth = 0.0;
	const index_t heavy_index = env.add_cell_type(heavy);

	for (index_t i = 0; i < 2 * per_type; ++i)
	{
		add_point_cell(env, i < per_type ? 0 : heavy_index, 10.0 * static_cast<real_t>(i));
	}

	openmp_solver solver;
	solver.initialize(env);
	solver.solve(env, steps);

	std::array<real_t, 2> mean_squared {};
	const auto positions = positions_of(env);
	for (index_t i = 0; i < 2 * per_type; ++i)
	{
		const real_t dx = positions[3 * i] - 10.0 * static_cast<real_t>(i);
		mean_squared.at(i < per_type ? 0 : 1) +=
			(dx * dx + positions[3 * i + 1] * positions[3 * i + 1] + positions[3 * i + 2] * positions[3 * i + 2])
			/ static_cast<real_t>(per_type);
	}

	// 2 * dims * D * time with D = thermal_energy / drag; 3000 elements give about 1.5% standard error
	const real_t time = timestep * static_cast<real_t>(steps);
	EXPECT_NEAR(mean_squared[0], 6.0 * thermal_energy / 1.0 * time, 0.08 * 6.0 * thermal_energy * time);
	EXPECT_NEAR(mean_squared[1], 6.0 * thermal_energy / 4.0 * time, 0.08 * 6.0 * thermal_energy / 4.0 * time);
}

// ---------------------------------------------------------------------------- walls and pressure

TEST(CellTypesTest, WallReachesHalfTheEquilibriumDistanceOfTheElementsType)
{
	environment env({ 3, { -10, -10, -10 }, { 10, 10, 10 }, { 1, 1, 1 } }, 0.001);
	env.params.domain_walls = true;
	const index_t large = env.add_cell_type(large_type());

	// both 0.8 from the wall at x = -10: beyond the reach of the small element (0.5), within that of
	// the large one (1.0)
	add_point_cell(env, 0, -9.2);
	add_point_cell(env, large, -9.2);
	env.agents->get_agent_at(1)->position()[1] = 5.0; // far from the first

	const auto forces = forces_of(env);
	EXPECT_DOUBLE_EQ(forces[0], 0.0);
	EXPECT_DOUBLE_EQ(forces[3], wall_force(wall_potential(env.cell_types[large], env.params), -9.2, -10.0, 10.0));
	EXPECT_GT(forces[3], 0.0);
}

TEST(CellTypesTest, PressureUsesTheVolumeAndTheWellsOfTheCellsType)
{
	// A cell of the large type with elements at -2 and +2, and outside it one element of the default
	// type 1.2 beyond each of them. The mixed equilibrium distance is 1.5, so they repel.
	const real_t gap = 1.2;
	environment env = make_environment(3);
	const index_t large = env.add_cell_type(large_type());

	env.create_cell(large);
	env.create_cell(0);
	env.create_agent(0)->position()[0] = -2.0;
	env.create_agent(0)->position()[0] = 2.0;
	env.create_agent(1)->position()[0] = -2.0 - gap;
	env.create_agent(1)->position()[0] = 2.0 + gap;

	cell_solver cells;
	cells.update(env);
	pressure_solver pressures;
	pressures.update(env);

	const interaction_table table(env.cell_types, env.params);
	const real_t coefficient = force_coefficient(table.inter(large, 0), gap * gap);
	ASSERT_GT(coefficient, 0.0);

	// each member is 2 from the centre and is pushed inwards with coefficient * gap; the cell has two
	// elements of the volume 2^3 / sqrt(2) each
	const real_t volume = 2.0 * 8.0 / std::numbers::sqrt2;
	const real_t expected = 2.0 * (2.0 * coefficient * gap) / (3.0 * volume);

	EXPECT_NEAR(env.cells->get_agent_at(0)->pressure(), expected, 1e-12 * expected);
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
