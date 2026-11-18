#include <array>
#include <cmath>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <string>

#include <gtest/gtest.h>

#include "micromechanics/environment_builder.h"
#include "micromechanics/solver_registry.h"

// The builder puts an environment together and refuses what cannot work, among it a timestep beyond
// the stability limit of forward Euler. These tests check the limit itself, against the formula and
// against what a cell actually does on either side of it.

namespace physicore::mechanics::micromechanics::tests {

namespace {
environment_builder builder_with_domain(index_t dims = 3)
{
	environment_builder builder;
	builder.resize(dims, { -50, -50, -50 }, { 50, 50, 50 }, { 10, 10, 10 });
	return builder;
}

// Largest distance of a sub-agent of the first cell from the centre of that cell
real_t reach_of_first_cell(environment& env)
{
	const auto centre = env.cells->get_agent_at(0)->position();
	real_t furthest = 0.0;
	for (index_t i = 0; i < env.agents->size(); ++i)
	{
		const auto position = env.agents->get_agent_at(i)->position();
		real_t squared = 0.0;
		for (std::size_t d = 0; d < position.size(); ++d)
		{
			squared += (position[d] - centre[d]) * (position[d] - centre[d]);
		}
		furthest = std::max(furthest, std::sqrt(squared));
	}
	return furthest;
}
} // namespace

TEST(StableTimestepTest, FollowsDragOverStiffnessWithTheFactorOfTheDimension)
{
	const simulation_parameters params; // steepness 2
	cell_type type;						// well depth 1, equilibrium distance 1, drag 1

	// pair stiffness 8 * 1 * 2^2 / 1^2 = 32
	EXPECT_DOUBLE_EQ(largest_stable_timestep(type, params, 1), 1.0 / (2.0 * 32.0));
	EXPECT_DOUBLE_EQ(largest_stable_timestep(type, params, 2), 1.0 / (3.0 * 32.0));
	EXPECT_DOUBLE_EQ(largest_stable_timestep(type, params, 3), 1.0 / (4.0 * 32.0));

	// twice the drag allows twice the step; four times the well depth a quarter; twice the size four times
	const real_t reference = largest_stable_timestep(type, params, 3);
	type.drag_coefficient = 2.0;
	EXPECT_DOUBLE_EQ(largest_stable_timestep(type, params, 3), 2.0 * reference);
	type.drag_coefficient = 1.0;
	type.intra_well_depth = 4.0;
	EXPECT_DOUBLE_EQ(largest_stable_timestep(type, params, 3), 0.25 * reference);
	type.intra_well_depth = 1.0;
	type.equilibrium_distance = 2.0;
	EXPECT_DOUBLE_EQ(largest_stable_timestep(type, params, 3), 4.0 * reference);
}

TEST(StableTimestepTest, TheDeeperOfTheTwoWellsCounts)
{
	const simulation_parameters params;
	cell_type type;
	const real_t reference = largest_stable_timestep(type, params, 3);

	type.inter_well_depth = 0.9; // still shallower than the intra-cellular well
	EXPECT_DOUBLE_EQ(largest_stable_timestep(type, params, 3), reference);

	type.inter_well_depth = 2.0;
	EXPECT_DOUBLE_EQ(largest_stable_timestep(type, params, 3), 0.5 * reference);
}

class StableTimestepDimsTest : public ::testing::TestWithParam<index_t>
{};

INSTANTIATE_TEST_SUITE_P(Dimensions, StableTimestepDimsTest, ::testing::Values(2, 3));

TEST_P(StableTimestepDimsTest, CellHoldsBelowTheLimitAndFliesApartAboveIt)
{
	const index_t dims = GetParam();
	const cell_type type;
	const real_t limit = largest_stable_timestep(type, simulation_parameters {}, dims);
	const std::array<real_t, 3> centre = { 0.0, 0.0, 0.0 };
	const index_t elements = dims == 2 ? 61 : 87; // four rings of a disc, five shells of a ball

	const auto reach_after = [&](real_t timestep) {
		// by hand and not with the builder, which would refuse the larger step
		environment env({ dims, { -50, -50, -50 }, { 50, 50, 50 }, { 10, 10, 10 } }, timestep);
		env.create_cell(std::span<const real_t>(centre.data(), dims), elements);
		env.solver = solver_registry::instance().get("openmp_solver");
		env.solver->initialize(env);
		const real_t at_start = reach_of_first_cell(env);
		try
		{
			env.solver->solve(env, 3000);
		}
		catch (const std::runtime_error&)
		{
			// the solver stops when a position is no longer finite
			return std::array<real_t, 2> { at_start, std::numeric_limits<real_t>::infinity() };
		}
		return std::array<real_t, 2> { at_start, reach_of_first_cell(env) };
	};

	// well below the limit the cell relaxes a little and stays what it is
	const auto stable = reach_after(0.5 * limit);
	EXPECT_LT(stable[1], 1.1 * stable[0]);
	EXPECT_GT(stable[1], 0.7 * stable[0]);

	// above it the elements oscillate with growing amplitude until they are out of each other's reach:
	// the cell ends scattered over a region much larger than itself
	const auto unstable = reach_after(2.0 * limit);
	EXPECT_GT(unstable[1], 1.5 * unstable[0]);
}

TEST(EnvironmentBuilderTest, BuildsAnEnvironmentThatIsReadyToRun)
{
	environment_builder builder = builder_with_domain();
	builder.set_time_step(0.001);
	simulation_parameters params;
	params.thermal_energy = 0.02;
	params.random_seed = 5;
	builder.set_parameters(params);
	builder.set_domain_walls(true);

	cell_type large;
	large.name = "large";
	large.equilibrium_distance = 2.0;
	large.elements_per_cell = 13;
	EXPECT_EQ(builder.add_cell_type(cell_type {}), index_t { 0 });
	EXPECT_EQ(builder.add_cell_type(large), index_t { 1 });

	auto env = builder.build();

	ASSERT_NE(env, nullptr);
	EXPECT_DOUBLE_EQ(env->mechanics_timestep, 0.001);
	EXPECT_EQ(env->mesh.dims, index_t { 3 });
	EXPECT_EQ(env->mesh.bounding_box_maxs[0], 50);
	EXPECT_TRUE(env->params.domain_walls);
	EXPECT_DOUBLE_EQ(env->params.thermal_energy, 0.02);
	ASSERT_EQ(env->cell_types.size(), 2U);
	EXPECT_EQ(env->cell_types[1].name, "large");
	ASSERT_NE(env->solver, nullptr);
	EXPECT_EQ(env->serializer, nullptr);

	const std::array<real_t, 3> centre = { 0.0, 0.0, 0.0 };
	env->create_cell_of_type(1, centre);
	EXPECT_NO_THROW(env->run_single_timestep());
	EXPECT_EQ(env->cells->get_agent_at(0)->agent_count(), index_t { 13 });
}

TEST(EnvironmentBuilderTest, WithoutCellTypesTheEnvironmentHasTheDefaultOne)
{
	environment_builder builder = builder_with_domain();
	builder.set_time_step(0.001);

	auto env = builder.build();

	ASSERT_EQ(env->cell_types.size(), 1U);
	EXPECT_DOUBLE_EQ(env->cell_types[0].equilibrium_distance, 1.0);
}

TEST(EnvironmentBuilderTest, RefusesToBuildWithoutDomainOrTimestep)
{
	environment_builder no_domain;
	no_domain.set_time_step(0.001);
	EXPECT_THROW(no_domain.build(), std::invalid_argument);

	environment_builder no_timestep = builder_with_domain();
	EXPECT_THROW(no_timestep.build(), std::invalid_argument);
}

TEST(EnvironmentBuilderTest, RefusesATimestepBeyondTheStabilityLimitAndSaysWhichTypeAndHowLarge)
{
	// the default type allows 1 / 128 = 0.0078125 in three dimensions
	environment_builder fine = builder_with_domain();
	fine.set_time_step(0.007);
	EXPECT_NO_THROW(fine.build());

	environment_builder builder = builder_with_domain();
	builder.set_time_step(0.007);
	builder.add_cell_type(cell_type {});
	cell_type stiff;
	stiff.name = "stiff";
	stiff.intra_well_depth = 4.0; // allows a quarter of that
	builder.add_cell_type(stiff);

	try
	{
		builder.build();
		FAIL() << "the timestep should have been refused";
	}
	catch (const std::invalid_argument& e)
	{
		const std::string message = e.what();
		EXPECT_NE(message.find("stiff"), std::string::npos) << message;
		EXPECT_NE(message.find("0.00195"), std::string::npos) << message;
	}

	// the same timestep is fine on a line, where the limit is twice as high... for the default type only
	environment_builder line = builder_with_domain(1);
	line.set_time_step(0.012);
	EXPECT_NO_THROW(line.build());
}

TEST(EnvironmentBuilderTest, RefusesInvalidCellTypesAndUnknownSolvers)
{
	environment_builder broken = builder_with_domain();
	broken.set_time_step(0.001);
	cell_type no_drag;
	no_drag.drag_coefficient = 0.0;
	broken.add_cell_type(no_drag);
	EXPECT_THROW(broken.build(), std::invalid_argument);

	environment_builder unknown = builder_with_domain();
	unknown.set_time_step(0.001);
	unknown.set_solver("no_such_solver");
	EXPECT_THROW(unknown.build(), std::invalid_argument);
}

TEST(EnvironmentBuilderTest, VtkOutputIsWrittenWhereAsked)
{
	const std::filesystem::path folder = std::filesystem::temp_directory_path() / "physicore_sem_builder_vtk";
	std::filesystem::remove_all(folder);

	environment_builder builder = builder_with_domain();
	builder.set_time_step(0.001);
	builder.set_vtk_output(folder.string());
	auto env = builder.build();
	const std::array<real_t, 3> centre = { 0.0, 0.0, 0.0 };
	env->create_cell(centre, 5);

	ASSERT_NE(env->serializer, nullptr);
	env->serialize_state(0.0);
	EXPECT_TRUE(std::filesystem::exists(folder / "micromechanics_agents.pvd"));

	std::filesystem::remove_all(folder);
}

} // namespace physicore::mechanics::micromechanics::tests
