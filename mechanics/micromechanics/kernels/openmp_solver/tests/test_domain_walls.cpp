#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <span>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>
#include <micromechanics/domain_wall.h>

#ifdef _OPENMP
	#include <omp.h>
#endif

#include "force_solver.h"
#include "openmp_solver.h"
#include "parallel_threshold.h"
#include "test_helpers.h"

// With domain_walls on, the faces of the domain push elements back in. An element at the distance h
// from a wall feels the repulsion of its own mirror image at the distance 2h, through the wall
// potential, for h below half the equilibrium distance. These tests check that force, that it keeps
// elements inside, and what it does to cells.

namespace physicore::mechanics::micromechanics::kernels::openmp_solver::tests {

namespace {
// A domain from -half_size to +half_size along every axis. The walls are switched on by the caller:
// an environment cannot be moved, so it cannot be changed here before it is returned.
environment make_box(index_t dims, real_t timestep = 0.001, sindex_t half_size = 10)
{
	return { { dims, { -half_size, -half_size, -half_size }, { half_size, half_size, half_size }, { 1, 1, 1 } },
			 timestep };
}

// Net force on every sub-agent, agents_count * dims values
std::vector<real_t> forces_of(environment& env)
{
	force_solver forces;
	forces.compute(env);
	return { forces.forces().begin(), forces.forces().end() };
}

// Total potential energy with the walls of a box of the given half size
real_t total_energy_with_walls(environment& env, real_t half_size)
{
	real_t energy = total_energy(env);
	for (index_t i = 0; i < env.agents->size(); ++i)
	{
		for (const real_t coordinate : env.agents->get_agent_at(i)->position())
		{
			energy += wall_energy(env.params.wall_potential(), coordinate, -half_size, half_size);
		}
	}
	return energy;
}
} // namespace

// ------------------------------------------------------------------------------------ the force

TEST(DomainWallTest, ForceIsTheRepulsionOfTheMirrorImage)
{
	const simulation_parameters params;
	const pair_potential wall = params.wall_potential();

	// 0.3 from the lower wall: the image is 0.6 away, closer than the equilibrium distance of 1
	const real_t expected = force_coefficient(wall, 0.6 * 0.6) * 0.6;
	ASSERT_GT(expected, 0.0);

	EXPECT_DOUBLE_EQ(wall_force(wall, 0.3, 0.0, 20.0), expected);
	// the upper wall pushes the other way with the same strength (20 - 0.3 is not exactly 0.3 from 20)
	EXPECT_NEAR(wall_force(wall, 20.0 - 0.3, 0.0, 20.0), -expected, 1e-9);
}

TEST(DomainWallTest, ForceActsOnlyWithinHalfAnEquilibriumDistanceAndNeverAttracts)
{
	const simulation_parameters params;
	const pair_potential wall = params.wall_potential();

	// touching: a sphere of radius r_eq / 2 just reaches the wall
	EXPECT_DOUBLE_EQ(wall_force(wall, -10.0 + 0.5, -10.0, 10.0), 0.0);
	// further in, where two elements at that distance would attract each other
	for (const real_t distance : { 0.55, 0.7, 1.0, 5.0 })
	{
		EXPECT_DOUBLE_EQ(wall_force(wall, -10.0 + distance, -10.0, 10.0), 0.0) << "distance " << distance;
	}

	// always away from the wall on the way in
	for (int k = 1; k < 50; ++k)
	{
		const real_t distance = 0.01 * k;
		EXPECT_GT(wall_force(wall, -10.0 + distance, -10.0, 10.0), 0.0) << "distance " << distance;
		EXPECT_LT(wall_force(wall, 10.0 - distance, -10.0, 10.0), 0.0) << "distance " << distance;
	}
}

TEST(DomainWallTest, ForceIsMinusTheDerivativeOfTheWallEnergy)
{
	const simulation_parameters params;
	const pair_potential wall = params.wall_potential();
	const real_t step = 1e-6;

	for (const real_t distance : { 0.05, 0.2, 0.35, 0.45 })
	{
		const real_t x = -10.0 + distance;
		const real_t derivative =
			(wall_energy(wall, x + step, -10.0, 10.0) - wall_energy(wall, x - step, -10.0, 10.0)) / (2.0 * step);
		EXPECT_NEAR(wall_force(wall, x, -10.0, 10.0), -derivative, 1e-6) << "distance " << distance;
	}

	// the energy is zero where the wall does not act and rises towards the wall
	EXPECT_DOUBLE_EQ(wall_energy(wall, 0.0, -10.0, 10.0), 0.0);
	EXPECT_NEAR(wall_energy(wall, -10.0 + 0.5, -10.0, 10.0), 0.0, 1e-15);
	EXPECT_GT(wall_energy(wall, -10.0 + 0.2, -10.0, 10.0), wall_energy(wall, -10.0 + 0.4, -10.0, 10.0));
}

TEST(DomainWallTest, ElementOnOrBeyondAWallFeelsNoForceFromIt)
{
	// the repulsion of this potential vanishes when two elements coincide; the solver puts such
	// elements back instead (see the tests of the hard stop)
	const simulation_parameters params;
	const pair_potential wall = params.wall_potential();

	EXPECT_DOUBLE_EQ(wall_force(wall, -10.0, -10.0, 10.0), 0.0);
	EXPECT_DOUBLE_EQ(wall_force(wall, -12.0, -10.0, 10.0), 0.0);
	EXPECT_DOUBLE_EQ(wall_force(wall, 10.5, -10.0, 10.0), 0.0);
}

class DomainWallDimsTest : public ::testing::TestWithParam<index_t>
{};

INSTANTIATE_TEST_SUITE_P(Dimensions, DomainWallDimsTest, ::testing::Values(1, 2, 3));

TEST_P(DomainWallDimsTest, EveryFaceOfTheDomainPushesAlongItsNormalOnly)
{
	const index_t dims = GetParam();
	const real_t expected = force_coefficient(simulation_parameters().wall_potential(), 0.6 * 0.6) * 0.6;

	for (index_t axis = 0; axis < dims; ++axis)
	{
		for (const real_t side : { -1.0, 1.0 })
		{
			environment env = make_box(dims);
			env.params.domain_walls = true;
			auto* element = env.agents->create();
			element->position()[axis] = side * (10.0 - 0.3);

			const auto forces = forces_of(env);
			for (index_t d = 0; d < dims; ++d)
			{
				EXPECT_NEAR(forces[d], d == axis ? -side * expected : 0.0, 1e-9)
					<< "wall on axis " << axis << " side " << side << ", force along " << d;
			}
		}
	}
}

TEST(DomainWallTest, TwoDimensionalDomainHasNoWallAlongZ)
{
	// the mesh has an extent along z, but a 2D simulation does not use that axis
	environment env = make_box(2);
	env.params.domain_walls = true;
	auto* element = env.agents->create();
	element->position()[0] = 3.0;
	element->position()[1] = -4.0;

	const auto forces = forces_of(env);
	ASSERT_EQ(forces.size(), 2U);
	EXPECT_DOUBLE_EQ(forces[0], 0.0);
	EXPECT_DOUBLE_EQ(forces[1], 0.0);
}

TEST(DomainWallTest, WallsAreOffByDefaultAndThenChangeNothing)
{
	const auto run = [](bool mention_walls) {
		// a tiny domain: everything is outside it
		environment env({ 3, { -1, -1, -1 }, { 1, 1, 1 }, { 1, 1, 1 } }, 0.002);
		if (mention_walls)
		{
			env.params.wall_well_depth = 7.0; // has no effect while the walls are off
		}
		EXPECT_FALSE(env.params.domain_walls);
		env.params.diffusion_coefficient = 0.01;
		add_jittered_lattice_cell(env, 1.0, 0.1, 5.0);

		openmp_solver solver;
		solver.initialize(env);
		solver.solve(env, 100);
		EXPECT_EQ(solver.wall_stops(), index_t { 0 });
		return positions_of(env);
	};

	EXPECT_EQ(run(false), run(true));
}

// ------------------------------------------------------------------------------ keeping elements in

TEST(DomainWallTest, CrawlingElementStopsWhereTheWallBalancesItsTraction)
{
	const real_t speed = 0.8;

	environment env = make_box(3, 0.0005);
	env.params.domain_walls = true;
	auto* cell = env.create_cell();
	cell->is_motile() = 1;
	cell->migration_speed() = speed;
	cell->persistence_time() = 1.0;
	cell->migration_bias() = 1.0;
	cell->migration_bias_direction()[0] = 1.0; // towards the wall at x = +10
	auto* element = add_agent(env, 0, { 8.0, 1.0, -2.0 });

	openmp_solver solver;
	solver.initialize(env);
	solver.solve(env, 20000);

	// At rest the wall force equals the traction drag_coefficient * speed. The wall force rises from
	// zero at the distance 0.5; the element stops at the first distance where it reaches the traction.
	const pair_potential wall = env.params.wall_potential();
	const real_t traction = env.params.drag_coefficient * speed;
	real_t far = 0.5;
	real_t near = 0.4;
	ASSERT_GT(-wall_force(wall, 10.0 - near, -10.0, 10.0), traction);
	for (int k = 0; k < 60; ++k)
	{
		const real_t middle = 0.5 * (far + near);
		(-wall_force(wall, 10.0 - middle, -10.0, 10.0) > traction ? near : far) = middle;
	}

	EXPECT_NEAR(10.0 - element->position()[0], 0.5 * (far + near), 1e-6);
	EXPECT_NEAR(element->velocity()[0], 0.0, 1e-9);
	EXPECT_DOUBLE_EQ(element->position()[1], 1.0);
	EXPECT_DOUBLE_EQ(element->position()[2], -2.0);
	EXPECT_EQ(solver.wall_stops(), index_t { 0 });
}

TEST_P(DomainWallDimsTest, NoElementLeavesTheDomainUnderStrongNoise)
{
	const index_t dims = GetParam();
	const sindex_t half_size = 3;

	// thermal energy D * drag = 2 is four times the well depth of the wall: kicks do carry elements
	// past the soft wall, and the hard stop has to catch them
	environment env = make_box(dims, 0.01, half_size);
	env.params.domain_walls = true;
	env.params.diffusion_coefficient = 2.0;
	env.params.inter_well_depth = 0.0;
	env.params.random_seed = 4;
	for (index_t i = 0; i < 200; ++i)
	{
		env.agents->create();
	}

	openmp_solver solver;
	solver.initialize(env);

	for (int block = 0; block < 20; ++block)
	{
		solver.solve(env, 100);
		for (const real_t coordinate : positions_of(env))
		{
			ASSERT_GE(coordinate, -static_cast<real_t>(half_size));
			ASSERT_LE(coordinate, static_cast<real_t>(half_size));
		}
	}

	// the elements have spread over the whole box, so they did reach the walls
	real_t furthest = 0.0;
	for (const real_t coordinate : positions_of(env))
	{
		furthest = std::max(furthest, std::abs(coordinate));
	}
	EXPECT_GT(furthest, static_cast<real_t>(half_size) - 0.5);
	EXPECT_GT(solver.wall_stops(), index_t { 0 });
}

TEST(DomainWallTest, GentleNoiseNeverNeedsTheHardStop)
{
	// thermal energy 0.02, a twenty-fifth of the well depth of the wall: the soft wall alone holds
	environment env = make_box(3, 0.002, 3);
	env.params.domain_walls = true;
	env.params.diffusion_coefficient = 0.02;
	env.params.inter_well_depth = 0.0;
	for (index_t i = 0; i < 100; ++i)
	{
		// start next to the walls, where it matters
		auto position = env.agents->create()->position();
		position[0] = 2.4;
		position[1] = -2.4;
		position[2] = 0.05 * static_cast<real_t>(i) - 2.4;
	}

	openmp_solver solver;
	solver.initialize(env);
	solver.solve(env, 5000);

	EXPECT_EQ(solver.wall_stops(), index_t { 0 });
	for (const real_t coordinate : positions_of(env))
	{
		EXPECT_LT(std::abs(coordinate), 3.0);
	}
}

TEST(DomainWallTest, ImmovableCellIsNotPutBack)
{
	environment env = make_box(3);
	env.params.domain_walls = true;
	env.create_cell();
	auto* element = add_agent(env, 0, { 0.0, 0.0, 0.0 });

	openmp_solver solver;
	solver.initialize(env);

	// moved outside by hand after the check at initialisation, and fixed there
	env.cells->get_agent_at(0)->is_movable() = 0;
	element->position()[0] = 25.0;
	solver.solve(env, 10);

	EXPECT_DOUBLE_EQ(element->position()[0], 25.0);
	EXPECT_EQ(solver.wall_stops(), index_t { 0 });
}

TEST(DomainWallTest, InitialisingWithAnElementOutsideIsAnError)
{
	environment env = make_box(3);
	env.params.domain_walls = true;
	env.agents->create()->position()[1] = 10.0; // on the wall counts as inside
	openmp_solver solver;
	EXPECT_NO_THROW(solver.initialize(env));

	env.agents->create()->position()[2] = -10.5;
	EXPECT_THROW(solver.initialize(env), std::invalid_argument);

	// without walls the size of the domain plays no role
	env.params.domain_walls = false;
	EXPECT_NO_THROW(solver.initialize(env));
}

// ----------------------------------------------------------------------------------------- cells

TEST(DomainWallTest, CellRelaxingAgainstAWallNeverGainsEnergy)
{
	// a cell placed so that some of its elements start inside the reach of the wall at x = -4
	const real_t half_size = 4.0;
	environment env = make_box(3, 0.0005, 4);
	env.params.domain_walls = true;
	add_jittered_lattice_cell(env, 0.95, 0.1, -3.8);

	openmp_solver solver;
	solver.initialize(env);

	real_t energy = total_energy_with_walls(env, half_size);
	const real_t initial_energy = energy;
	for (int block = 0; block < 200; ++block)
	{
		solver.solve(env, 10);
		const real_t now = total_energy_with_walls(env, half_size);
		ASSERT_LE(now, energy + 1e-12) << "block " << block;
		energy = now;
	}

	EXPECT_LT(energy, initial_energy - 0.01);
	EXPECT_EQ(solver.wall_stops(), index_t { 0 });
	// the wall pushed the cell inwards
	for (index_t i = 0; i < env.agents->size(); ++i)
	{
		EXPECT_GT(env.agents->get_agent_at(i)->position()[0], -half_size + 0.3);
	}
}

TEST(DomainWallTest, CellCrawlingIntoAWallIsCompressedByItsTraction)
{
	const real_t speed = 0.2;
	const index_t elements = 13;

	environment env = make_box(3, 0.002, 6);
	env.params.domain_walls = true;
	const std::array<real_t, 3> centre = { 3.0, 0.0, 0.0 };
	auto* cell = env.create_cell(centre, elements);
	cell->is_motile() = 1;
	cell->migration_speed() = speed;
	cell->persistence_time() = 1.0;
	cell->migration_bias() = 1.0;
	cell->migration_bias_direction()[0] = 1.0;

	openmp_solver solver;
	solver.initialize(env);
	EXPECT_DOUBLE_EQ(cell->pressure(), 0.0); // not touching anything yet

	solver.solve(env, 15000);

	// stalled against the wall at x = 6
	EXPECT_NEAR(cell->velocity()[0], 0.0, 0.01 * speed);
	EXPECT_GT(cell->position()[0], 4.0);
	EXPECT_LT(cell->position()[0], 6.0);
	EXPECT_EQ(solver.wall_stops(), index_t { 0 });

	// The wall holds the cell with a force equal to its traction N * drag * speed, applied at the wall,
	// which is the distance from the centre of the cell to the wall away. A little less in fact: the
	// force acts on elements up to half an equilibrium distance inside the wall.
	const real_t traction = static_cast<real_t>(elements) * env.params.drag_coefficient * speed;
	const real_t lever_arm = 6.0 - cell->position()[0];
	const real_t volume = static_cast<real_t>(elements) / std::numbers::sqrt2;
	const real_t estimate = traction * lever_arm / (3.0 * volume);

	EXPECT_GT(cell->pressure(), 0.0);
	EXPECT_NEAR(cell->pressure(), estimate, 0.4 * estimate);
}

TEST(DomainWallTest, GrowingCellInACornerStaysInsideTheDomain)
{
	environment env = make_box(2, 0.002, 3);
	env.params.domain_walls = true;
	const std::array<real_t, 2> centre = { 2.0, 2.0 };
	env.create_cell(centre, 7);

	for (int k = 0; k < 40; ++k)
	{
		const auto position = env.grow_cell(0)->position();
		ASSERT_LE(std::abs(position[0]), 3.0) << "growth " << k;
		ASSERT_LE(std::abs(position[1]), 3.0) << "growth " << k;
	}

	// without walls the same cell does grow beyond that box
	environment open = make_box(2, 0.002, 3);
	open.create_cell(centre, 7);
	real_t furthest = 0.0;
	for (int k = 0; k < 40; ++k)
	{
		const auto position = open.grow_cell(0)->position();
		furthest = std::max({ furthest, std::abs(position[0]), std::abs(position[1]) });
	}
	EXPECT_GT(furthest, 3.0);
}

// -------------------------------------------------------------------------------- reproducibility

#ifdef _OPENMP
TEST(DomainWallTest, RunWithWallsDoesNotDependOnTheNumberOfThreads)
{
	const auto run = [](int threads) {
		omp_set_num_threads(threads);

		// a small box and strong noise, so that the walls and the hard stop are both at work
		environment env = make_box(3, 0.01, 2);
		env.params.domain_walls = true;
		env.params.diffusion_coefficient = 2.0;
		env.params.inter_well_depth = 0.0;
		for (index_t i = 0; i < 4 * parallel_threshold; ++i)
		{
			env.agents->create();
		}

		openmp_solver solver;
		solver.initialize(env);
		solver.solve(env, 60);

		std::vector<real_t> result = positions_of(env);
		result.push_back(static_cast<real_t>(solver.wall_stops()));
		return result;
	};

	const int original_threads = omp_get_max_threads();
	const auto single_threaded = run(1);
	const auto multi_threaded = run(4);
	omp_set_num_threads(original_threads);

	EXPECT_EQ(single_threaded, multi_threaded);
	EXPECT_GT(single_threaded.back(), 0.0);
}
#endif

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver::tests
