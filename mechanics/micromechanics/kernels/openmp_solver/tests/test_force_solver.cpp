#include <array>
#include <cmath>
#include <random>
#include <vector>

#include <gtest/gtest.h>

#ifdef _OPENMP
	#include <omp.h>
#endif

#include "force_solver.h"
#include "test_helpers.h"

namespace physicore::mechanics::micromechanics::kernels::openmp_solver::tests {

namespace {
// Deliberately not unit values, so a wrong power of a parameter cannot cancel out.
// Cutoff distance = 2.5 * 2 = 5.
void set_test_parameters(environment& env)
{
	env.params.equilibrium_distance = 2.0;
	env.params.steepness = 1.5;
	env.params.intra_well_depth = 3.0;
	env.params.inter_well_depth = 0.75;
	env.params.cutoff_factor = 2.5;
}

std::vector<real_t> computed_forces(environment& env)
{
	force_solver solver;
	solver.compute(env);
	const auto forces = solver.forces();
	return { forces.begin(), forces.end() };
}
} // namespace

TEST(ForceSolverTest, TwoElementsOfOneCellUseTheIntraPotential)
{
	environment env = make_environment(3);
	set_test_parameters(env);
	env.create_cell();
	add_agent(env, 0, { 0.0, 0.0, 0.0 });
	add_agent(env, 0, { 1.5, 0.0, 0.0 });

	const auto forces = computed_forces(env);

	// closer than the equilibrium distance: repulsion pushes each element away from the other
	const real_t c = force_coefficient(env.params.intra_potential(), 2.25);
	ASSERT_GT(c, 0.0);
	ASSERT_EQ(forces.size(), 6U);
	EXPECT_DOUBLE_EQ(forces[0], -1.5 * c);
	EXPECT_DOUBLE_EQ(forces[3], 1.5 * c);
	for (const std::size_t k : { 1U, 2U, 4U, 5U })
	{
		EXPECT_DOUBLE_EQ(forces[k], 0.0);
	}
}

TEST(ForceSolverTest, ElementsOfDifferentCellsUseTheInterPotential)
{
	environment env = make_environment(3);
	set_test_parameters(env);
	env.create_cell();
	env.create_cell();
	add_agent(env, 0, { 0.0, 0.0, 0.0 });
	add_agent(env, 1, { 0.0, 3.0, 0.0 });

	const auto forces = computed_forces(env);

	// further apart than the equilibrium distance: adhesion pulls the elements together
	const real_t c = force_coefficient(env.params.inter_potential(), 9.0);
	ASSERT_LT(c, 0.0);
	EXPECT_DOUBLE_EQ(forces[1], -3.0 * c);
	EXPECT_DOUBLE_EQ(forces[4], 3.0 * c);
	EXPECT_GT(forces[1], 0.0);

	// the inter-cellular well is four times shallower than the intra-cellular one
	EXPECT_DOUBLE_EQ(4.0 * c, force_coefficient(env.params.intra_potential(), 9.0));
}

TEST(ForceSolverTest, ElementsWithoutACellInteractThroughTheInterPotential)
{
	environment env = make_environment(3);
	set_test_parameters(env);
	env.agents->create()->position()[0] = 0.0;
	env.agents->create()->position()[0] = 1.0;

	const auto forces = computed_forces(env);

	EXPECT_DOUBLE_EQ(forces[3], force_coefficient(env.params.inter_potential(), 1.0));
}

TEST(ForceSolverTest, NoForceFromTheCutoffOn)
{
	environment env = make_environment(3);
	set_test_parameters(env);
	env.create_cell();
	add_agent(env, 0, { 0.0, 0.0, 0.0 });
	add_agent(env, 0, { 5.0, 0.0, 0.0 });
	add_agent(env, 0, { 0.0, -40.0, 0.0 });

	for (const real_t force : computed_forces(env))
	{
		EXPECT_DOUBLE_EQ(force, 0.0);
	}
}

TEST(ForceSolverTest, HandComputedTriangle)
{
	// Right triangle with both legs at the equilibrium distance: the legs carry no force and only the
	// hypotenuse pair (distance^2 = 8) interacts.
	environment env = make_environment(2);
	set_test_parameters(env);
	env.create_cell();
	add_agent(env, 0, { 0.0, 0.0 });
	add_agent(env, 0, { 2.0, 0.0 });
	add_agent(env, 0, { 0.0, 2.0 });

	const auto forces = computed_forces(env);

	const real_t c = force_coefficient(env.params.intra_potential(), 8.0);
	ASSERT_LT(c, 0.0);
	ASSERT_EQ(forces.size(), 6U);
	EXPECT_DOUBLE_EQ(forces[0], 0.0);
	EXPECT_DOUBLE_EQ(forces[1], 0.0);
	EXPECT_DOUBLE_EQ(forces[2], 2.0 * c);
	EXPECT_DOUBLE_EQ(forces[3], -2.0 * c);
	EXPECT_DOUBLE_EQ(forces[4], -2.0 * c);
	EXPECT_DOUBLE_EQ(forces[5], 2.0 * c);
}

TEST(ForceSolverTest, PlanarConfigurationGivesTheSameForcesIn2DAnd3D)
{
	const std::vector<std::array<real_t, 2>> points = { { 0.0, 0.0 }, { 1.7, 0.3 }, { -0.4, 2.2 }, { 2.5, 2.9 } };

	environment planar = make_environment(2);
	environment spatial = make_environment(3);
	for (environment* env : { &planar, &spatial })
	{
		set_test_parameters(*env);
		env->create_cell();
		env->create_cell();
	}
	for (std::size_t i = 0; i < points.size(); ++i)
	{
		add_agent(planar, i % 2, { points[i][0], points[i][1] });
		add_agent(spatial, i % 2, { points[i][0], points[i][1], 0.0 });
	}

	const auto planar_forces = computed_forces(planar);
	const auto spatial_forces = computed_forces(spatial);

	for (std::size_t i = 0; i < points.size(); ++i)
	{
		EXPECT_DOUBLE_EQ(spatial_forces[3 * i], planar_forces[2 * i]);
		EXPECT_DOUBLE_EQ(spatial_forces[3 * i + 1], planar_forces[2 * i + 1]);
		EXPECT_DOUBLE_EQ(spatial_forces[3 * i + 2], 0.0);
		EXPECT_NE(planar_forces[2 * i], 0.0);
	}
}

TEST(ForceSolverTest, NetForceAndNetTorqueVanish)
{
	// Pair forces are equal, opposite and act along the line joining the two elements, so an isolated
	// system feels neither a net force nor a net torque, whatever the mix of cells.
	environment env = make_environment(3);
	for (index_t i = 0; i < 3; ++i)
	{
		env.create_cell();
	}

	std::mt19937 generator(7);
	std::uniform_real_distribution<real_t> coordinate(0.0, 4.0);
	for (index_t i = 0; i < 60; ++i)
	{
		add_agent(env, i % 3, { coordinate(generator), coordinate(generator), coordinate(generator) });
	}

	const auto forces = computed_forces(env);
	const auto positions = positions_of(env);

	std::array<real_t, 3> net_force {};
	std::array<real_t, 3> net_torque {};
	real_t largest_force = 0.0;
	for (std::size_t i = 0; i < 60; ++i)
	{
		const real_t* x = &positions[3 * i];
		const real_t* f = &forces[3 * i];

		for (std::size_t d = 0; d < 3; ++d)
		{
			net_force.at(d) += f[d];
			largest_force = std::max(largest_force, std::abs(f[d]));
		}
		net_torque[0] += x[1] * f[2] - x[2] * f[1];
		net_torque[1] += x[2] * f[0] - x[0] * f[2];
		net_torque[2] += x[0] * f[1] - x[1] * f[0];
	}

	ASSERT_GT(largest_force, 1.0);
	for (std::size_t d = 0; d < 3; ++d)
	{
		EXPECT_NEAR(net_force.at(d), 0.0, 1e-10 * largest_force);
		EXPECT_NEAR(net_torque.at(d), 0.0, 1e-10 * largest_force);
	}
}

#ifdef _OPENMP
TEST(ForceSolverTest, ResultDoesNotDependOnTheNumberOfThreads)
{
	environment env = make_environment(3);
	add_jittered_lattice_cell(env, 1.0, 0.2);

	const int original_threads = omp_get_max_threads();

	omp_set_num_threads(1);
	const auto single_threaded = computed_forces(env);
	omp_set_num_threads(4);
	const auto multi_threaded = computed_forces(env);

	omp_set_num_threads(original_threads);

	// bitwise identical, not merely close
	EXPECT_EQ(single_threaded, multi_threaded);
}
#endif

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver::tests
