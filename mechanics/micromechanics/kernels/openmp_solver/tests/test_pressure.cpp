#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <span>
#include <vector>

#include <gtest/gtest.h>

#ifdef _OPENMP
	#include <omp.h>
#endif

#include "cell_solver.h"
#include "openmp_solver.h"
#include "parallel_threshold.h"
#include "pressure_solver.h"
#include "test_helpers.h"

// The pressure of a cell is the virial of the forces that sub-agents outside the cell exert on it,
//
//     pressure = -1 / (dims * volume) * sum (x_i - X_cell) . F_ij
//
// which returns p for a cell whose surface is pressed with a uniform pressure p. These tests check the
// value against cases computed by hand, its sign, its symmetries and the grid against the brute-force sum.

namespace physicore::mechanics::micromechanics::kernels::openmp_solver::tests {

namespace {
// volume of one sub-agent in a close packing at the default equilibrium distance of 1
real_t element_volume(index_t dims)
{
	if (dims == 1)
		return 1.0;
	return dims == 2 ? 0.5 * std::numbers::sqrt3 : 1.0 / std::numbers::sqrt2;
}

void update_pressures(environment& env)
{
	cell_solver cells;
	cells.update(env);

	pressure_solver pressures;
	pressures.update(env);
}

std::vector<real_t> pressures_of(environment& env)
{
	std::vector<real_t> pressures;
	pressures.reserve(env.cells->size());
	for (index_t c = 0; c < env.cells->size(); ++c)
	{
		pressures.push_back(env.cells->get_agent_at(c)->pressure());
	}
	return pressures;
}

// A cell of two sub-agents at -1 and +1 on the x axis, and outside it one sub-agent of another cell
// at the distance `gap` beyond each of them. The outer sub-agents are 2 + gap from the far member of the
// cell, which is beyond the cutoff of 2.5 for the gaps used here, so each acts on one member only.
void add_squeezed_pair(environment& env, real_t gap)
{
	env.create_cell();
	env.create_cell();

	env.create_agent(0)->position()[0] = -1.0;
	env.create_agent(0)->position()[0] = 1.0;
	env.create_agent(1)->position()[0] = -1.0 - gap;
	env.create_agent(1)->position()[0] = 1.0 + gap;
}
} // namespace

class PressureTest : public ::testing::TestWithParam<index_t>
{};

INSTANTIATE_TEST_SUITE_P(Dimensions, PressureTest, ::testing::Values(1, 2, 3));

// cells made by create_cell, which packs them in two or three dimensions
class PackedCellPressureTest : public ::testing::TestWithParam<index_t>
{};

INSTANTIATE_TEST_SUITE_P(Dimensions, PackedCellPressureTest, ::testing::Values(2, 3));

TEST_P(PressureTest, MatchesTheValueComputedByHand)
{
	const index_t dims = GetParam();
	const real_t gap = 0.8;

	environment env = make_environment(dims);
	add_squeezed_pair(env, gap);
	update_pressures(env);

	// The member at +1 is pushed by its neighbour at 1 + gap with the force c * (x_i - x_j) = -c * gap
	// along x, c being the force coefficient of the inter-cellular potential at that distance. It sits
	// at +1 from the centre of its cell, so its term of the virial is -c * gap. The member at -1 gives
	// the same. The cell has two sub-agents.
	const real_t coefficient = force_coefficient(inter_potential(env.cell_types[0], env.params), gap * gap);
	ASSERT_GT(coefficient, 0.0); // closer than the equilibrium distance: repulsion
	const real_t expected = 2.0 * coefficient * gap / (static_cast<real_t>(dims) * 2.0 * element_volume(dims));

	EXPECT_NEAR(env.cells->get_agent_at(0)->pressure(), expected, 1e-12 * expected);
	EXPECT_GT(env.cells->get_agent_at(0)->pressure(), 0.0);
}

TEST_P(PressureTest, IsNegativeWhenTheNeighboursPullOutwards)
{
	const index_t dims = GetParam();
	const real_t gap = 1.2; // further than the equilibrium distance: attraction

	environment env = make_environment(dims);
	add_squeezed_pair(env, gap);
	update_pressures(env);

	const real_t coefficient = force_coefficient(inter_potential(env.cell_types[0], env.params), gap * gap);
	ASSERT_LT(coefficient, 0.0);
	const real_t expected = 2.0 * coefficient * gap / (static_cast<real_t>(dims) * 2.0 * element_volume(dims));

	EXPECT_NEAR(env.cells->get_agent_at(0)->pressure(), expected, -1e-12 * expected);
	EXPECT_LT(env.cells->get_agent_at(0)->pressure(), 0.0);
}

TEST_P(PressureTest, IsZeroAtTheEquilibriumDistanceAndBeyondTheCutoff)
{
	const index_t dims = GetParam();

	environment touching = make_environment(dims);
	add_squeezed_pair(touching, 1.0);
	update_pressures(touching);
	EXPECT_NEAR(touching.cells->get_agent_at(0)->pressure(), 0.0, 1e-14);

	environment apart = make_environment(dims);
	add_squeezed_pair(apart, 2.6);
	update_pressures(apart);
	EXPECT_DOUBLE_EQ(apart.cells->get_agent_at(0)->pressure(), 0.0);
}

TEST_P(PackedCellPressureTest, IsolatedCellHasNoPressureWhateverItsInternalForces)
{
	const index_t dims = GetParam();
	environment env = make_environment(dims);
	const std::array<real_t, 3> centre = { 3.0, -2.0, 1.0 };
	// packed at the equilibrium distance with the default steepness the cell is not at rest: its
	// second neighbours attract. Those forces are inside the cell and must not count.
	env.create_cell(std::span<const real_t>(centre.data(), dims), 37);
	update_pressures(env);

	EXPECT_DOUBLE_EQ(env.cells->get_agent_at(0)->pressure(), 0.0);
}

TEST(PressureTest, SubAgentsWithoutACellPressOnACellButHaveNoPressureThemselves)
{
	// as the squeezed pair, but the outer sub-agents belong to no cell
	const real_t gap = 0.8;
	environment env = make_environment(3);
	env.create_cell();
	env.create_agent(0)->position()[0] = -1.0;
	env.create_agent(0)->position()[0] = 1.0;
	env.agents->create()->position()[0] = -1.0 - gap;
	env.agents->create()->position()[0] = 1.0 + gap;

	update_pressures(env);

	environment reference = make_environment(3);
	add_squeezed_pair(reference, gap);
	update_pressures(reference);

	EXPECT_DOUBLE_EQ(env.cells->get_agent_at(0)->pressure(), reference.cells->get_agent_at(0)->pressure());
}

TEST(PressureTest, CellOfOneSubAgentAndEmptyCellHaveNoPressure)
{
	environment env = make_environment(3);
	env.create_cell();
	env.create_cell();
	env.create_cell();
	add_agent(env, 0, { 0.0, 0.0, 0.0 });
	// pressed hard by a neighbour, but a single sub-agent sits on the centre of its cell: no lever arm
	add_agent(env, 1, { 0.7, 0.0, 0.0 });

	update_pressures(env);

	EXPECT_DOUBLE_EQ(env.cells->get_agent_at(0)->pressure(), 0.0);
	EXPECT_DOUBLE_EQ(env.cells->get_agent_at(1)->pressure(), 0.0);
	EXPECT_DOUBLE_EQ(env.cells->get_agent_at(2)->pressure(), 0.0);
}

TEST(PressureTest, PressureGoesBackToZeroWhenTheNeighbourLeaves)
{
	environment env = make_environment(3);
	add_squeezed_pair(env, 0.8);
	update_pressures(env);
	ASSERT_GT(env.cells->get_agent_at(0)->pressure(), 0.0);

	env.remove_cell(1);
	update_pressures(env);

	EXPECT_DOUBLE_EQ(env.cells->get_agent_at(0)->pressure(), 0.0);
}

TEST_P(PackedCellPressureTest, TwoMirroredCellsHaveTheSamePressure)
{
	const index_t dims = GetParam();
	// The packing of create_cell is symmetric under inversion, so two cells placed at -x and +x are
	// mirror images of each other through the origin.
	for (const real_t half_distance : { 1.0, 1.4 })
	{
		environment env = make_environment(dims);
		for (const real_t side : { -1.0, 1.0 })
		{
			const std::array<real_t, 3> centre = { side * half_distance, 0.0, 0.0 };
			env.create_cell(std::span<const real_t>(centre.data(), dims), 19);
		}
		update_pressures(env);

		const auto pressures = pressures_of(env);
		EXPECT_NE(pressures[0], 0.0);
		EXPECT_NEAR(pressures[0], pressures[1], 1e-12 * std::abs(pressures[0])) << "half distance " << half_distance;
	}
}

TEST(PressureTest, OverlappingCellsAreCompressedAndCellsHeldAtADistanceAreUnderTension)
{
	const auto pressure_at = [](real_t distance) {
		environment env = make_environment(3);
		for (const real_t side : { -0.5, 0.5 })
		{
			const std::array<real_t, 3> centre = { side * distance, 0.0, 0.0 };
			env.create_cell(centre, 13);
		}
		update_pressures(env);
		return env.cells->get_agent_at(0)->pressure();
	};

	// The twelve outer sub-agents of a cell sit one equilibrium distance from its centre, the nearest
	// ones 0.71 along x. At a centre distance of 2.0 the facing sub-agents of the two cells are 0.59
	// apart and repel; at 2.8 they are 1.39 apart and attract.
	EXPECT_GT(pressure_at(2.0), 0.0);
	EXPECT_LT(pressure_at(2.8), 0.0);
	// beyond reach: the facing sub-agents are 3.2 apart
	EXPECT_DOUBLE_EQ(pressure_at(4.6), 0.0);
}

TEST(PressureTest, DoesNotDependOnWhereTheSystemIs)
{
	const auto run = [](real_t shift) {
		environment env = make_environment(3);
		add_jittered_lattice_cell(env, 1.0, 0.1, shift, 42);
		add_jittered_lattice_cell(env, 1.0, 0.1, shift + 3.1, 43);
		add_jittered_lattice_cell(env, 1.0, 0.1, shift + 6.0, 44);
		update_pressures(env);
		return pressures_of(env);
	};

	const auto here = run(0.0);
	const auto there = run(57.25);

	for (std::size_t c = 0; c < here.size(); ++c)
	{
		EXPECT_NE(here[c], 0.0);
		EXPECT_NEAR(here[c], there[c], 1e-9 * std::abs(here[c])) << "cell " << c;
	}
}

TEST(PressureTest, GridAgreesWithTheBruteForceSum)
{
	environment env = make_environment(3);
	for (int c = 0; c < 6; ++c)
	{
		add_jittered_lattice_cell(env, 1.0, 0.15, 2.9 * c, 100 + c);
	}
	// and a few loose sub-agents between them
	for (int k = 0; k < 5; ++k)
	{
		auto position = env.agents->create()->position();
		position[0] = 2.9 * k + 2.45;
		position[1] = 1.0;
		position[2] = 1.0;
	}

	cell_solver cells;
	cells.update(env);

	pressure_solver pressures;
	pressures.update(env);
	const auto with_grid = pressures_of(env);

	pressures.update_brute_force(env);
	const auto brute_force = pressures_of(env);

	for (std::size_t c = 0; c < with_grid.size(); ++c)
	{
		EXPECT_NE(brute_force[c], 0.0);
		EXPECT_NEAR(with_grid[c], brute_force[c], 1e-10 * std::abs(brute_force[c])) << "cell " << c;
	}
}

TEST(PressureTest, SolverKeepsThePressureUpToDate)
{
	environment env = make_environment(3, 0.002);
	for (const real_t side : { -1.0, 1.0 })
	{
		const std::array<real_t, 3> centre = { 1.0 * side, 0.0, 0.0 };
		env.create_cell(centre, 13);
	}
	auto* left = env.cells->get_agent_at(0);

	openmp_solver solver;
	solver.initialize(env);

	// the cells start overlapping: compressed
	const real_t at_start = left->pressure();
	EXPECT_GT(at_start, 0.0);

	// they push each other apart and settle in adhesive contact, where push and pull nearly cancel
	solver.solve(env, 5000);

	EXPECT_LT(std::abs(left->pressure()), 0.05 * at_start);
	EXPECT_NEAR(left->velocity()[0], 0.0, 1e-3);
}

TEST(PressureTest, CellsCrawlingAgainstEachOtherAreCompressedByTheirTraction)
{
	// Two motile cells driven head-on along a line stall against each other (see the motility tests).
	// Each is then pressed by the other with its total traction N * drag * speed, acting on its last
	// sub-agent. The virial of that is the traction times the distance of that sub-agent from the
	// centre of its cell.
	const real_t speed = 0.2;
	const index_t elements = 5;

	environment env = make_environment(1, 0.002);
	for (const real_t side : { -1.0, 1.0 })
	{
		const std::array<real_t, 1> centre = { 5.0 * side };
		auto* cell = env.create_cell(centre, elements);
		cell->is_motile() = 1;
		cell->migration_speed() = speed;
		cell->persistence_time() = 1.0;
		cell->migration_bias() = 1.0;
		cell->migration_bias_direction()[0] = -side;
	}
	auto* left = env.cells->get_agent_at(0);
	auto* right = env.cells->get_agent_at(1);

	openmp_solver solver;
	solver.initialize(env);
	solver.solve(env, 20000);
	ASSERT_NEAR(left->velocity()[0], 0.0, 0.002);

	real_t rightmost_of_left = -1e30;
	for (index_t i = 0; i < env.agents->size(); ++i)
	{
		auto* element = env.agents->get_agent_at(i);
		if (element->cell_index() == 0)
			rightmost_of_left = std::max(rightmost_of_left, element->position()[0]);
	}

	const real_t traction = static_cast<real_t>(elements) * env.cell_types[0].drag_coefficient * speed;
	const real_t lever_arm = rightmost_of_left - left->position()[0];
	const real_t volume = static_cast<real_t>(elements) * element_volume(1);
	const real_t estimate = traction * lever_arm / (1.0 * volume);

	EXPECT_GT(left->pressure(), 0.0);
	// not exact: the second sub-agent from the end is also in reach of the other cell and is pulled by it
	EXPECT_NEAR(left->pressure(), estimate, 0.1 * estimate);
	EXPECT_NEAR(left->pressure(), right->pressure(), 1e-6 * left->pressure());
}

TEST(PressureTest, PressureUsesTheVolumeAndTheWellsOfTheCellsType)
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

TEST(PressureTest, TypesWithIdenticalValuesBehaveLikeOneType)
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

TEST(PressureTest, CellsWithoutAffinityStillCannotPassThroughEachOther)
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

#ifdef _OPENMP
TEST(PressureTest, DoesNotDependOnTheNumberOfThreads)
{
	const auto run = [](int threads) {
		omp_set_num_threads(threads);

		environment env = make_environment(3);
		// 12 cells of 27 sub-agents: more than the parallel threshold
		for (int c = 0; c < 12; ++c)
		{
			add_jittered_lattice_cell(env, 1.0, 0.15, 2.9 * c, 200 + c);
		}
		EXPECT_GT(env.agents->size(), parallel_threshold);
		update_pressures(env);
		return pressures_of(env);
	};

	const int original_threads = omp_get_max_threads();
	const auto single_threaded = run(1);
	const auto multi_threaded = run(4);
	omp_set_num_threads(original_threads);

	EXPECT_EQ(single_threaded, multi_threaded);
}
#endif

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver::tests
