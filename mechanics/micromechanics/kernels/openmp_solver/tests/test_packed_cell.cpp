#include <algorithm>
#include <array>
#include <cmath>

#include <gtest/gtest.h>
#include <micromechanics/cell_parameters.h>

#include "openmp_solver.h"
#include "test_helpers.h"

// A cell created with environment::create_cell(centre, elements) starts without overlapping elements.
// That matters for stiff potentials: elements that start on top of each other are thrown apart, and
// beyond the cutoff nothing pulls them back.

namespace physicore::mechanics::micromechanics::kernels::openmp_solver::tests {

namespace {
constexpr cell_parameters paper_cell { .radius = 10.0, .stiffness = 5000.0, .damping = 5000.0 };

struct cell_extent
{
	real_t furthest = 0.0;
	real_t radius_from_gyration = 0.0;
};

cell_extent extent_of(environment& env)
{
	const auto centre = centre_of_mass(env);
	const index_t count = env.agents->size();

	cell_extent extent;
	real_t squared_gyration_radius = 0.0;
	for (index_t i = 0; i < count; ++i)
	{
		const auto position = env.agents->get_agent_at(i)->position();
		real_t squared = 0.0;
		for (std::size_t d = 0; d < 3; ++d)
		{
			squared += std::pow(position[d] - centre.at(d), 2);
		}
		extent.furthest = std::max(extent.furthest, std::sqrt(squared));
		squared_gyration_radius += squared / static_cast<real_t>(count);
	}
	extent.radius_from_gyration = std::sqrt(5.0 / 3.0 * squared_gyration_radius);
	return extent;
}
} // namespace

class PackedCellTest : public ::testing::TestWithParam<real_t>
{};

INSTANTIATE_TEST_SUITE_P(Steepness, PackedCellTest, ::testing::Values(2.0, 4.0));

TEST_P(PackedCellTest, PackedCellStaysIntactAndKeepsTheRequestedSize)
{
	simulation_parameters shape;
	shape.steepness = GetParam();
	const simulation_parameters params = scale_to_elements(paper_cell, 256, shape);

	const real_t pair_stiffness = 8.0 * params.intra_well_depth * params.steepness * params.steepness
								  / (params.equilibrium_distance * params.equilibrium_distance);
	const real_t timestep = 0.05 * params.drag_coefficient / pair_stiffness;

	environment env = make_environment(3, timestep);
	env.params = params;
	const std::array<real_t, 3> centre = { 3.0, -2.0, 7.0 };
	env.create_cell(centre, 256);

	const cell_extent before = extent_of(env);
	const auto initial_centre = centre_of_mass(env);

	// two seconds, twice the relaxation time of the cell
	openmp_solver solver;
	solver.solve(env, static_cast<index_t>(2.0 / timestep));

	const cell_extent after = extent_of(env);

	// nobody is thrown out: the cell only settles, it does not expand
	EXPECT_LT(after.furthest, 1.05 * before.furthest);
	// The size stays within a tenth of the requested radius. Measured: 9.23 for steepness 2, which
	// compresses itself, and 10.00 for steepness 4, which keeps the close packing it started in.
	EXPECT_NEAR(after.radius_from_gyration, paper_cell.radius, 0.1 * paper_cell.radius);
	// and the cell does not drift
	const auto final_centre = centre_of_mass(env);
	for (std::size_t d = 0; d < 3; ++d)
	{
		EXPECT_NEAR(final_centre.at(d), initial_centre.at(d), 1e-9);
	}
}

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver::tests
