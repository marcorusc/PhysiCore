#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#ifdef _OPENMP
	#include <omp.h>
#endif

#include "openmp_solver.h"
#include "parallel_threshold.h"
#include "test_helpers.h"

// Motility is a persistent random walk of the cell (the model of PhysiCell): every sub-agent of a
// motile cell gets the extra velocity migration_speed * motility_direction, and the direction is
// redrawn with probability timestep / persistence_time per step. These tests check the motion of one
// cell, the statistics of the walk, reproducibility and what happens when two cells meet.

namespace physicore::mechanics::micromechanics::kernels::openmp_solver::tests {

namespace {
constexpr real_t separation = 10.0;

// Motile cells of one sub-agent each that do not interact: the inter-cellular potential is switched
// off. They start on a line along x. All get the same motility parameters and no direction yet.
void add_walkers(environment& env, index_t count, real_t speed, real_t persistence_time, real_t bias = 0.0)
{
	env.params.inter_well_depth = 0.0;
	for (index_t c = 0; c < count; ++c)
	{
		auto* walker = env.create_cell();
		walker->is_motile() = 1;
		walker->migration_speed() = speed;
		walker->persistence_time() = persistence_time;
		walker->migration_bias() = bias;

		env.create_agent(c)->position()[0] = separation * static_cast<real_t>(c);
	}
}

// How far every walker has moved from where add_walkers() put it, cells_count * dims values
std::vector<real_t> displacements_of(environment& env)
{
	const index_t dims = env.mesh.dims;

	std::vector<real_t> displacements = positions_of(env);
	for (index_t i = 0; i < env.agents->size(); ++i)
	{
		displacements[i * dims] -= separation * static_cast<real_t>(i);
	}
	return displacements;
}

// Copy of all motility directions, cells_count * dims values
std::vector<real_t> directions_of(environment& env)
{
	std::vector<real_t> directions;
	for (index_t c = 0; c < env.cells->size(); ++c)
	{
		const auto direction = env.cells->get_agent_at(c)->motility_direction();
		directions.insert(directions.end(), direction.begin(), direction.end());
	}
	return directions;
}

// Mean squared displacement of an unbiased walker after `steps` steps, exact for the discrete walk.
// The displacement is speed * timestep * (d_0 + ... + d_{steps-1}). Two directions m steps apart are
// the same if the cell did not turn in between, which has probability q^m with q = 1 - timestep /
// persistence_time, and are independent with zero mean otherwise, so <d_j . d_k> = q^|j-k|.
real_t discrete_mean_squared_displacement(real_t speed, real_t persistence_time, real_t timestep, index_t steps)
{
	const real_t q = 1.0 - timestep / persistence_time;

	auto sum = static_cast<real_t>(steps);
	real_t power = 1.0;
	for (index_t m = 1; m < steps; ++m)
	{
		power *= q;
		sum += 2.0 * static_cast<real_t>(steps - m) * power;
	}
	return speed * speed * timestep * timestep * sum;
}

// The same in continuous time (Fuerth's formula): ballistic at short times, diffusive at long times
real_t continuous_mean_squared_displacement(real_t speed, real_t persistence_time, real_t time)
{
	return 2.0 * speed * speed * persistence_time
		   * (time - persistence_time * (1.0 - std::exp(-time / persistence_time)));
}
} // namespace

// ------------------------------------------------------------------------------------- one cell

class MotileCellTest : public ::testing::TestWithParam<index_t>
{};

INSTANTIATE_TEST_SUITE_P(Dimensions, MotileCellTest, ::testing::Values(2, 3));

TEST_P(MotileCellTest, IsolatedCellMovesAtTheMigrationSpeedWithoutChangingShape)
{
	const index_t dims = GetParam();
	const real_t timestep = 0.002;
	const index_t steps = 500;
	const real_t speed = 0.7;
	// a unit vector that is not along an axis
	const std::array<real_t, 3> heading =
		dims == 2 ? std::array<real_t, 3> { 0.6, -0.8, 0.0 } : std::array<real_t, 3> { 0.6, 0.0, -0.8 };

	// The same cell twice, once motile and once not. The cell is not at rest (it relaxes during the
	// run), so the motile run is compared with the other one, not with the start.
	const auto run = [&](bool motile) {
		environment env = make_environment(dims, timestep);
		const std::array<real_t, 3> centre = { 1.0, 2.0, 3.0 };
		auto* cell = env.create_cell(std::span<const real_t>(centre.data(), dims), 19);
		cell->is_motile() = motile ? 1 : 0;
		cell->migration_speed() = speed;
		cell->persistence_time() = 0.05; // turns often, but with a bias of 1 always onto the same heading
		cell->migration_bias() = 1.0;
		std::ranges::copy_n(heading.begin(), static_cast<std::ptrdiff_t>(dims),
							cell->migration_bias_direction().begin());

		openmp_solver solver;
		solver.initialize(env);
		solver.solve(env, steps);

		std::vector<real_t> result = positions_of(env);
		const auto velocity = cell->velocity();
		result.insert(result.end(), velocity.begin(), velocity.end());
		return result;
	};

	const auto resting = run(false);
	const auto crawling = run(true);

	const real_t time = timestep * static_cast<real_t>(steps);
	const index_t elements = 19;

	// every sub-agent is displaced by exactly speed * time * heading: the cell moved as a rigid body
	for (index_t i = 0; i < elements; ++i)
	{
		for (index_t d = 0; d < dims; ++d)
		{
			EXPECT_NEAR(crawling[i * dims + d] - resting[i * dims + d], speed * time * heading.at(d), 1e-10)
				<< "sub-agent " << i << " axis " << d;
		}
	}

	// the velocity of the cell is the migration velocity; the forces inside a cell cancel in their sum
	for (index_t d = 0; d < dims; ++d)
	{
		EXPECT_NEAR(crawling[elements * dims + d], speed * heading.at(d), 1e-10) << "axis " << d;
		EXPECT_NEAR(resting[elements * dims + d], 0.0, 1e-10) << "axis " << d;
	}
}

TEST(MotilityTest, CellThatIsNotMotileIgnoresItsMotilityParameters)
{
	environment env = make_environment(3, 0.01);
	auto* cell = env.create_cell();
	cell->migration_speed() = 5.0;
	cell->persistence_time() = 1.0;
	cell->migration_bias() = 1.0;
	cell->migration_bias_direction()[0] = 1.0;
	cell->motility_direction()[1] = 1.0;
	auto* element = add_agent(env, 0, { 1.0, 2.0, 3.0 });

	openmp_solver solver;
	solver.solve(env, 100);

	EXPECT_DOUBLE_EQ(element->position()[0], 1.0);
	EXPECT_DOUBLE_EQ(element->position()[1], 2.0);
	EXPECT_DOUBLE_EQ(element->position()[2], 3.0);
	// and its direction is left alone
	EXPECT_DOUBLE_EQ(cell->motility_direction()[0], 0.0);
	EXPECT_DOUBLE_EQ(cell->motility_direction()[1], 1.0);
}

TEST(MotilityTest, ImmovableCellDoesNotCrawl)
{
	environment env = make_environment(3, 0.01);
	add_walkers(env, 1, 5.0, 1.0);
	env.cells->get_agent_at(0)->is_movable() = 0;

	openmp_solver solver;
	solver.solve(env, 100);

	for (const real_t x : displacements_of(env))
	{
		EXPECT_DOUBLE_EQ(x, 0.0);
	}
	EXPECT_DOUBLE_EQ(env.agents->get_agent_at(0)->velocity()[0], 0.0);
}

TEST(MotilityTest, SubAgentWithoutCellDoesNotCrawl)
{
	environment env = make_environment(3, 0.01);
	add_walkers(env, 1, 5.0, 1.0);
	auto* loose = env.agents->create();
	loose->position()[1] = 50.0;

	openmp_solver solver;
	solver.solve(env, 10);

	EXPECT_DOUBLE_EQ(loose->position()[0], 0.0);
	EXPECT_DOUBLE_EQ(loose->position()[1], 50.0);
	EXPECT_DOUBLE_EQ(loose->position()[2], 0.0);
}

// ------------------------------------------------------------------------------- the random walk

class PersistentWalkTest : public ::testing::TestWithParam<index_t>
{};

INSTANTIATE_TEST_SUITE_P(Dimensions, PersistentWalkTest, ::testing::Values(1, 2, 3));

TEST_P(PersistentWalkTest, MeanSquaredDisplacementFollowsThePersistentRandomWalk)
{
	const index_t dims = GetParam();
	const index_t walkers = 4000;
	const real_t speed = 1.5;
	const real_t persistence_time = 0.5;
	const real_t timestep = 0.01;
	const index_t steps = 200; // four persistence times: from ballistic into diffusive motion

	environment env = make_environment(dims, timestep);
	env.params.random_seed = 11;
	add_walkers(env, walkers, speed, persistence_time);

	openmp_solver solver;

	// measured halfway and at the end, to see the shape of the curve and not just one value
	for (const index_t steps_so_far : { steps / 2, steps })
	{
		solver.solve(env, steps / 2);

		real_t measured = 0.0;
		for (const real_t x : displacements_of(env))
		{
			measured += x * x / static_cast<real_t>(walkers);
		}

		const real_t expected = discrete_mean_squared_displacement(speed, persistence_time, timestep, steps_so_far);

		// 4000 independent walkers: the standard error of the estimate is a few percent at most
		EXPECT_NEAR(measured, expected, 0.1 * expected) << "after " << steps_so_far << " steps";

		// the walk in discrete steps is the continuous one up to a correction of order
		// timestep / persistence_time = 2%
		const real_t time = timestep * static_cast<real_t>(steps_so_far);
		EXPECT_NEAR(expected, continuous_mean_squared_displacement(speed, persistence_time, time), 0.03 * expected);
	}

	// an unbiased walk goes nowhere on average: the mean displacement along each axis is zero within
	// five standard errors
	const auto displacements = displacements_of(env);
	const real_t per_axis =
		discrete_mean_squared_displacement(speed, persistence_time, timestep, steps) / static_cast<real_t>(dims);
	const real_t standard_error = std::sqrt(per_axis / static_cast<real_t>(walkers));
	for (index_t d = 0; d < dims; ++d)
	{
		real_t mean = 0.0;
		for (index_t c = 0; c < walkers; ++c)
		{
			mean += displacements[c * dims + d] / static_cast<real_t>(walkers);
		}
		EXPECT_NEAR(mean, 0.0, 5.0 * standard_error) << "axis " << d;
	}
}

TEST(MotilityTest, StraightRunsLastThePersistenceTimeOnAverage)
{
	const index_t walkers = 2000;
	const index_t steps = 500;
	const real_t timestep = 0.01;
	const real_t persistence_time = 0.5;

	environment env = make_environment(3, timestep);
	add_walkers(env, walkers, 0.0, persistence_time);

	openmp_solver solver;
	solver.solve(env, 1); // every walker draws its first direction

	// an unbiased walker in 3D never draws the same direction twice, so every turn shows
	index_t turns = 0;
	auto before = directions_of(env);
	for (index_t s = 0; s < steps; ++s)
	{
		solver.solve(env, 1);
		const auto after = directions_of(env);
		for (index_t c = 0; c < walkers; ++c)
		{
			if (after[3 * c] != before[3 * c])
			{
				++turns;
			}
		}
		before = after;
	}

	// a walker turns with probability p = timestep / persistence_time = 0.02 per step. In a million
	// trials the measured fraction has a standard error of sqrt(p (1 - p) / 1e6) = 1.4e-4.
	const real_t turn_fraction = static_cast<real_t>(turns) / static_cast<real_t>(walkers * steps);
	EXPECT_NEAR(turn_fraction, timestep / persistence_time, 5.0 * 1.4e-4);

	// the mean number of steps between two turns is 1 / p, so the mean run lasts the persistence time
	EXPECT_NEAR(timestep / turn_fraction, persistence_time, 0.04 * persistence_time);
}

TEST_P(PersistentWalkTest, UnbiasedDirectionsAreUnitVectorsWithoutPreferredAxis)
{
	const index_t dims = GetParam();
	const index_t walkers = 20000;

	environment env = make_environment(dims, 0.01);
	add_walkers(env, walkers, 0.0, 1.0);

	openmp_solver solver;
	solver.solve(env, 1);

	const auto directions = directions_of(env);
	std::array<real_t, 3> mean {};
	std::array<real_t, 3> mean_square {};
	for (index_t c = 0; c < walkers; ++c)
	{
		real_t squared_length = 0.0;
		for (index_t d = 0; d < dims; ++d)
		{
			const real_t component = directions[c * dims + d];
			squared_length += component * component;
			mean.at(d) += component / static_cast<real_t>(walkers);
			mean_square.at(d) += component * component / static_cast<real_t>(walkers);
		}
		ASSERT_NEAR(squared_length, 1.0, 1e-12) << "cell " << c;
	}

	// A component of a uniform unit vector has mean 0 and variance 1 / dims, so the mean of 20000 has
	// a standard error of sqrt(1 / (dims * 20000)), at most 0.0071.
	for (index_t d = 0; d < dims; ++d)
	{
		EXPECT_NEAR(mean.at(d), 0.0, 5.0 * 0.0071) << "axis " << d;
		EXPECT_NEAR(mean_square.at(d), 1.0 / static_cast<real_t>(dims), 0.02) << "axis " << d;
	}
}

TEST(MotilityTest, FullBiasFollowsTheBiasDirectionAndPartialBiasLeansTowardsIt)
{
	const index_t walkers = 5000;
	const std::array<real_t, 3> bias_direction = { 0.0, 0.6, 0.8 };

	const auto mean_direction = [&](real_t bias) {
		environment env = make_environment(3, 0.01);
		add_walkers(env, walkers, 0.0, 1.0, bias);
		for (index_t c = 0; c < walkers; ++c)
		{
			std::ranges::copy(bias_direction, env.cells->get_agent_at(c)->migration_bias_direction().begin());
		}

		openmp_solver solver;
		solver.solve(env, 1);

		const auto directions = directions_of(env);
		std::array<real_t, 3> mean {};
		for (index_t c = 0; c < walkers; ++c)
		{
			for (std::size_t d = 0; d < 3; ++d)
			{
				mean.at(d) += directions[3 * c + d] / static_cast<real_t>(walkers);
			}
		}
		return mean;
	};

	// how far the mean direction points along the bias direction: 1 if all cells follow it exactly
	const auto along_bias = [&](const std::array<real_t, 3>& mean) {
		return mean[0] * bias_direction[0] + mean[1] * bias_direction[1] + mean[2] * bias_direction[2];
	};

	EXPECT_NEAR(along_bias(mean_direction(1.0)), 1.0, 1e-12);

	const real_t unbiased = along_bias(mean_direction(0.0));
	const real_t weak = along_bias(mean_direction(0.25));
	const real_t strong = along_bias(mean_direction(0.75));

	EXPECT_NEAR(unbiased, 0.0, 5.0 * std::sqrt(1.0 / (3.0 * static_cast<real_t>(walkers))));
	EXPECT_GT(weak, 0.1);
	EXPECT_GT(strong, weak + 0.1);
	EXPECT_LT(strong, 1.0);
}

TEST(MotilityTest, DirectionIsKeptWhenBiasAndBiasDirectionGiveNoDirection)
{
	// a bias of 1 towards a bias direction that was never set: the mix has no length
	environment env = make_environment(3, 0.01);
	add_walkers(env, 1, 1.0, 0.001, 1.0); // would turn in every step
	auto* cell = env.cells->get_agent_at(0);
	cell->motility_direction()[2] = 1.0;

	openmp_solver solver;
	solver.solve(env, 50);

	EXPECT_DOUBLE_EQ(cell->motility_direction()[0], 0.0);
	EXPECT_DOUBLE_EQ(cell->motility_direction()[1], 0.0);
	EXPECT_DOUBLE_EQ(cell->motility_direction()[2], 1.0);
	EXPECT_NEAR(env.agents->get_agent_at(0)->position()[2], 50 * 0.01 * 1.0, 1e-12);
}

// ----------------------------------------------------------------------------------- parameters

TEST(MotilityTest, MotileCellNeedsAPositivePersistenceTimeAndABiasBetweenZeroAndOne)
{
	environment env = make_environment(3, 0.01);
	add_walkers(env, 2, 1.0, 1.0);
	auto* cell = env.cells->get_agent_at(1);

	openmp_solver solver;
	EXPECT_NO_THROW(solver.initialize(env));

	cell->persistence_time() = 0.0;
	EXPECT_THROW(solver.initialize(env), std::invalid_argument);
	EXPECT_THROW(solver.solve(env, 1), std::invalid_argument);

	cell->persistence_time() = 1.0;
	cell->migration_bias() = 1.5;
	EXPECT_THROW(solver.initialize(env), std::invalid_argument);

	cell->migration_bias() = -0.1;
	EXPECT_THROW(solver.solve(env, 1), std::invalid_argument);

	// a cell that is not motile is not checked
	cell->is_motile() = 0;
	EXPECT_NO_THROW(solver.solve(env, 1));
}

// -------------------------------------------------------------------------------- reproducibility

TEST(MotilityTest, SameSeedGivesTheSameWalkAndAnotherSeedADifferentOne)
{
	const auto run = [](std::uint64_t seed) {
		environment env = make_environment(3, 0.01);
		env.params.random_seed = seed;
		add_walkers(env, 50, 1.0, 0.1, 0.3);

		openmp_solver solver;
		solver.solve(env, 200);
		return positions_of(env);
	};

	const auto first = run(5);

	// bitwise identical, not merely close
	EXPECT_EQ(run(5), first);
	EXPECT_NE(run(6), first);
}

TEST(MotilityTest, MotilityDoesNotChangeTheThermalNoise)
{
	// Walkers with zero speed draw directions but do not crawl. Their thermal motion must be exactly
	// that of cells that are not motile: the two kinds of random numbers do not share anything.
	const auto run = [](bool motile) {
		environment env = make_environment(3, 0.01);
		env.params.diffusion_coefficient = 0.1;
		env.params.random_seed = 3;
		add_walkers(env, 50, 0.0, 0.1);
		for (index_t c = 0; c < 50; ++c)
		{
			env.cells->get_agent_at(c)->is_motile() = motile ? 1 : 0;
		}

		openmp_solver solver;
		solver.solve(env, 100);
		return positions_of(env);
	};

	EXPECT_EQ(run(true), run(false));
}

#ifdef _OPENMP
TEST(MotilityTest, WalkDoesNotDependOnTheNumberOfThreads)
{
	const auto run = [](int threads) {
		omp_set_num_threads(threads);

		environment env = make_environment(3, 0.01);
		add_walkers(env, 4 * parallel_threshold, 1.0, 0.1, 0.3);

		openmp_solver solver;
		solver.solve(env, 50);
		return positions_of(env);
	};

	const int original_threads = omp_get_max_threads();
	const auto single_threaded = run(1);
	const auto multi_threaded = run(4);
	omp_set_num_threads(original_threads);

	EXPECT_EQ(single_threaded, multi_threaded);
}
#endif

// ------------------------------------------------------------------------------------ two cells

TEST(MotilityTest, TwoCellsCrawlingHeadOnAlongALineStallAgainstEachOther)
{
	// Each sub-agent is pushed with the force drag_coefficient * migration_speed. The cells stop where
	// the repulsion between them balances that push.
	//
	// The cells are on a line (one dimension), where they cannot get past each other. In two or three
	// dimensions the head-on position is an unstable balance: the smallest sideways offset grows and
	// the cells slide past each other, so whether they stay jammed depends on rounding.
	const real_t speed = 0.2;
	const real_t timestep = 0.002;
	const index_t elements = 5;

	environment env = make_environment(1, timestep);
	for (const real_t side : { -1.0, 1.0 })
	{
		const std::array<real_t, 1> centre = { 5.0 * side };
		auto* cell = env.create_cell(centre, elements);
		cell->is_motile() = 1;
		cell->migration_speed() = speed;
		cell->persistence_time() = 1.0;
		cell->migration_bias() = 1.0;
		cell->migration_bias_direction()[0] = -side; // towards the other cell
	}
	auto* left = env.cells->get_agent_at(0);
	auto* right = env.cells->get_agent_at(1);

	openmp_solver solver;
	solver.initialize(env);
	const real_t centre_before = centre_of_mass(env)[0];
	const real_t separation_before = right->position()[0] - left->position()[0];

	// they first approach at the migration speed
	solver.solve(env, 50);
	EXPECT_NEAR(left->velocity()[0], speed, 0.05 * speed);
	EXPECT_NEAR(right->velocity()[0], -speed, 0.05 * speed);

	solver.solve(env, 20000);

	// stalled: both cells are at rest although they keep pushing
	EXPECT_NEAR(left->velocity()[0], 0.0, 0.01 * speed);
	EXPECT_NEAR(right->velocity()[0], 0.0, 0.01 * speed);

	// closer than at the start, and still in the same order: every sub-agent of the left cell is to
	// the left of every sub-agent of the right cell
	EXPECT_LT(right->position()[0] - left->position()[0], separation_before);
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
	// pressed together: closer than the equilibrium distance, but not through each other
	EXPECT_GT(leftmost_of_right - rightmost_of_left, 0.5);
	EXPECT_LT(leftmost_of_right - rightmost_of_left, 1.0);

	// the two pushes are equal and opposite and the forces between sub-agents cancel in pairs, so the
	// centre of all sub-agents does not move
	EXPECT_NEAR(centre_of_mass(env)[0], centre_before, 1e-9);
}

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver::tests
