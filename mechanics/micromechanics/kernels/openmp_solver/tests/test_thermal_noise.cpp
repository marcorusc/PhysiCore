#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

#include <gtest/gtest.h>

#ifdef _OPENMP
	#include <omp.h>
#endif

#include "openmp_solver.h"
#include "parallel_threshold.h"
#include "test_helpers.h"

// Thermal noise turns the equation of motion into the overdamped Langevin equation of the paper
// (Sandersius and Newman 2008, Eq. 1 and 2). These tests check its two physical signatures, free
// diffusion and thermal equilibrium, and that a noisy run is still reproducible.

namespace physicore::mechanics::micromechanics::kernels::openmp_solver::tests {

namespace {
// Elements that do not interact: without a cell they use the inter-cellular potential, which is
// switched off. They start on a line, far enough apart that the neighbour search has nothing to do.
void add_free_elements(environment& env, index_t count)
{
	env.cell_types[0].inter_well_depth = 0.0;
	for (index_t i = 0; i < count; ++i)
	{
		env.agents->create()->position()[0] = 10.0 * static_cast<real_t>(i);
	}
}

// How far every element has moved from where add_free_elements() put it, agents_count * dims values
std::vector<real_t> displacements_of(environment& env)
{
	const index_t dims = env.mesh.dims;

	std::vector<real_t> displacements = positions_of(env);
	for (index_t i = 0; i < env.agents->size(); ++i)
	{
		displacements[i * dims] -= 10.0 * static_cast<real_t>(i);
	}
	return displacements;
}
} // namespace

class FreeDiffusionTest : public ::testing::TestWithParam<index_t>
{};

INSTANTIATE_TEST_SUITE_P(Dimensions, FreeDiffusionTest, ::testing::Values(1, 2, 3));

TEST_P(FreeDiffusionTest, MeanSquaredDisplacementGrowsAsTwoDimsDTime)
{
	const index_t dims = GetParam();
	const index_t elements = 4000;
	const real_t diffusion = 0.16;
	const real_t timestep = 0.01;
	const index_t steps = 200;

	environment env = make_environment(dims, timestep);
	env.params.thermal_energy = diffusion;
	add_free_elements(env, elements);

	openmp_solver solver;

	// measured halfway and at the end, to see the growth with time and not just one value
	std::vector<real_t> mean_squared_displacement;
	for (int half = 0; half < 2; ++half)
	{
		solver.solve(env, steps / 2);

		real_t sum = 0.0;
		for (const real_t x : displacements_of(env))
		{
			sum += x * x;
		}
		mean_squared_displacement.push_back(sum / static_cast<real_t>(elements));
	}

	const real_t end_time = timestep * static_cast<real_t>(steps);
	const real_t expected = 2.0 * static_cast<real_t>(dims) * diffusion * end_time;

	// the estimate from 4000 elements has a relative standard error of sqrt(2 / (4000 * dims)), at most 2.3%
	EXPECT_NEAR(mean_squared_displacement[1], expected, 0.1 * expected);
	EXPECT_NEAR(mean_squared_displacement[0], 0.5 * expected, 0.1 * 0.5 * expected);

	// no drift: the mean displacement along each axis is zero within five standard errors
	const real_t standard_error = std::sqrt(2.0 * diffusion * end_time / static_cast<real_t>(elements));
	const auto displacements = displacements_of(env);
	for (index_t d = 0; d < dims; ++d)
	{
		real_t mean = 0.0;
		for (index_t i = 0; i < elements; ++i)
		{
			mean += displacements[i * dims + d] / static_cast<real_t>(elements);
		}
		EXPECT_NEAR(mean, 0.0, 5.0 * standard_error) << "axis " << d;
	}
}

TEST(ThermalNoiseTest, DisplacementScalesWithTheSquareRootOfTheTimestep)
{
	// One step of pure noise: variance of the displacement along one axis is 2 D dt.
	const index_t elements = 20000;
	const real_t diffusion = 0.5;

	for (const real_t timestep : { 0.001, 0.004 })
	{
		environment env = make_environment(3, timestep);
		env.params.thermal_energy = diffusion;
		add_free_elements(env, elements);

		openmp_solver solver;
		solver.solve(env, 1);

		real_t variance = 0.0;
		for (const real_t x : displacements_of(env))
		{
			variance += x * x / static_cast<real_t>(3 * elements);
		}
		EXPECT_NEAR(variance, 2.0 * diffusion * timestep, 0.03 * 2.0 * diffusion * timestep);
	}
}

TEST(ThermalNoiseTest, PairFluctuatesWithTheThermalEnergy)
{
	// In thermal equilibrium a spring of stiffness k stores the energy kT / 2 on average, so the
	// distance of a bound pair fluctuates with variance kT / k. For the overdamped Langevin equation
	// kT = D * drag_coefficient (Einstein relation).
	//
	// kT / k is the leading term. The measured variance is larger by three known corrections: the
	// finite timestep (factor 1 / (1 - k dt / drag), +2% here), the well being softer on its outer
	// side (+ kT * (36 g3^2 / k^3 - 12 g4 / k^2) with the cubic and quartic coefficients g3 = -48 and
	// g4 = 57.3 of the default potential) and the pair being free to rotate in three dimensions.
	// At kT = 1e-3 k they add up to about +10%, which is what a run at that temperature gives (+11%).
	// The test uses a temperature five times lower, where they add up to about +4% (measured: +3.6%).
	const index_t pairs = 400;
	const real_t stiffness = 32.0; // 8 * u0 * rho^2 / r_eq^2 with the default parameters
	const real_t thermal_energy = 2e-4 * stiffness;

	// relaxation time of the pair distance is drag / (2 k) = 1 / 64; the timestep resolves it 25 times
	const real_t timestep = 0.02 / stiffness;

	environment env = make_environment(3, timestep);
	env.params.thermal_energy = thermal_energy;
	env.params.random_seed = 2026;
	env.create_cell();

	// pairs far beyond the cutoff from each other; all in one cell, so they use the intra-cellular well
	for (index_t p = 0; p < pairs; ++p)
	{
		const real_t x = 10.0 * static_cast<real_t>(p);
		add_agent(env, 0, { x, 0.0, 0.0 });
		add_agent(env, 0, { x, 1.0, 0.0 });
	}

	openmp_solver solver;
	solver.solve(env, 500); // twenty relaxation times to forget the start

	real_t sum = 0.0;
	real_t sum_of_squares = 0.0;
	index_t samples = 0;
	for (int sample = 0; sample < 40; ++sample)
	{
		solver.solve(env, 50); // two relaxation times between samples
		for (index_t p = 0; p < pairs; ++p)
		{
			const real_t distance = std::sqrt(squared_distance(env, 2 * p, 2 * p + 1));
			sum += distance;
			sum_of_squares += distance * distance;
			++samples;
		}
	}

	const real_t mean = sum / static_cast<real_t>(samples);
	const real_t variance = sum_of_squares / static_cast<real_t>(samples) - mean * mean;

	EXPECT_NEAR(variance, thermal_energy / stiffness, 0.1 * thermal_energy / stiffness);
	// the pair stays bound at its equilibrium distance
	EXPECT_NEAR(mean, 1.0, 0.01);
}

TEST(ThermalNoiseTest, SameSeedGivesTheSameRunAndAnotherSeedADifferentOne)
{
	const auto run = [](std::uint64_t seed) {
		environment env = make_environment(3, 0.002);
		env.params.thermal_energy = 0.01;
		env.params.random_seed = seed;
		add_jittered_lattice_cell(env, 1.0, 0.1);

		openmp_solver solver;
		solver.solve(env, 200);
		return positions_of(env);
	};

	const auto first = run(5);

	// bitwise identical, not merely close
	EXPECT_EQ(run(5), first);
	EXPECT_NE(run(6), first);
}

TEST(ThermalNoiseTest, ZeroDiffusionIsTheDeterministicSolverWhateverTheSeed)
{
	const auto run = [](std::uint64_t seed) {
		environment env = make_environment(3, 0.002);
		env.params.random_seed = seed;
		add_jittered_lattice_cell(env, 1.0, 0.1);

		openmp_solver solver;
		solver.solve(env, 50);
		return positions_of(env);
	};

	EXPECT_EQ(run(1), run(987654321));
}

TEST(ThermalNoiseTest, ImmovableCellIsNotShakenByNoise)
{
	environment env = make_environment(3, 0.01);
	env.params.thermal_energy = 1.0;
	env.create_cell()->is_movable() = 0;
	auto* element = add_agent(env, 0, { 1.0, 2.0, 3.0 });

	openmp_solver solver;
	solver.solve(env, 100);

	EXPECT_DOUBLE_EQ(element->position()[0], 1.0);
	EXPECT_DOUBLE_EQ(element->position()[1], 2.0);
	EXPECT_DOUBLE_EQ(element->position()[2], 3.0);
}

TEST(ThermalNoiseTest, CellDiffusesMoreSlowlyThanItsElementsAndStaysIntact)
{
	// The forces inside a cell cancel in their sum, so the centre of a cell of N elements moves only
	// by the average of N independent noise terms: it diffuses with D / N. Many cells are followed
	// to measure that. Each is the origin and its twelve nearest neighbours, far from all others.
	const index_t cells_count = 200;
	const index_t elements_per_cell = 13;
	const real_t timestep = 0.002;
	const index_t steps = 200;

	environment env = make_environment(3, timestep);
	// thermal energy D * drag = 0.05 u0: strong enough to shake the cells, too weak to break bonds
	env.params.thermal_energy = 0.05;
	env.params.random_seed = 7;

	for (index_t c = 0; c < cells_count; ++c)
	{
		const std::array<real_t, 3> centre = { 10.0 * static_cast<real_t>(c), 0.0, 0.0 };
		env.create_cell(centre, elements_per_cell);
	}

	openmp_solver solver;
	solver.initialize(env);
	solver.solve(env, steps);

	real_t mean_squared_displacement = 0.0;
	for (index_t c = 0; c < cells_count; ++c)
	{
		auto* cell = env.cells->get_agent_at(c);
		ASSERT_EQ(cell->agent_count(), elements_per_cell);

		const std::array<real_t, 3> start = { 10.0 * static_cast<real_t>(c), 0.0, 0.0 };
		for (std::size_t d = 0; d < 3; ++d)
		{
			mean_squared_displacement +=
				std::pow(cell->position()[d] - start.at(d), 2) / static_cast<real_t>(cells_count);
		}
	}

	const real_t end_time = timestep * static_cast<real_t>(steps);
	const real_t expected = 2.0 * 3.0 * env.params.thermal_energy / static_cast<real_t>(elements_per_cell) * end_time;

	// 200 cells in 3D: relative standard error sqrt(2 / 600) = 6%. A longer run with 400 cells and
	// 500 steps measured 1.04 times the expected value.
	EXPECT_NEAR(mean_squared_displacement, expected, 0.2 * expected);

	// every element is still within reach of its cell: at creation the furthest is one equilibrium
	// distance from the centre
	for (index_t i = 0; i < env.agents->size(); ++i)
	{
		auto* element = env.agents->get_agent_at(i);
		auto* cell = env.cells->get_agent_at(element->cell_index());

		real_t squared = 0.0;
		for (std::size_t d = 0; d < 3; ++d)
		{
			squared += std::pow(element->position()[d] - cell->position()[d], 2);
		}
		ASSERT_LT(std::sqrt(squared), 1.6) << "element " << i;
	}
}

#ifdef _OPENMP
TEST(ThermalNoiseTest, NoisyRunDoesNotDependOnTheNumberOfThreads)
{
	const auto run = [](int threads) {
		omp_set_num_threads(threads);

		environment env = make_environment(3, 0.01);
		env.params.thermal_energy = 0.05;
		add_free_elements(env, 4 * parallel_threshold);

		openmp_solver solver;
		solver.solve(env, 20);
		return positions_of(env);
	};

	const int original_threads = omp_get_max_threads();
	const auto single_threaded = run(1);
	const auto multi_threaded = run(4);
	omp_set_num_threads(original_threads);

	EXPECT_EQ(single_threaded, multi_threaded);
}
#endif

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver::tests
