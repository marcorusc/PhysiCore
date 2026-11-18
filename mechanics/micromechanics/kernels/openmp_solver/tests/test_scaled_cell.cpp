#include <cmath>
#include <random>

#include <gtest/gtest.h>
#include <micromechanics/cell_parameters.h>

#include "openmp_solver.h"
#include "test_helpers.h"

// A cell is described by its radius, stiffness and damping, and scale_to_elements() turns those into
// element parameters for a chosen number of elements. These tests check the point of that scaling:
// the cell that emerges has the requested size, whatever the number of elements.

namespace physicore::mechanics::micromechanics::kernels::openmp_solver::tests {

namespace {
// The cell of Sandersius and Newman (2008), section 4.1, in micrometres, seconds and piconewtons
constexpr cell_parameters paper_cell { .radius = 10.0, .stiffness = 5000.0, .damping = 5000.0 };

// Radius of the relaxed cell: sqrt(5/3) times the radius of gyration, which is exact for a uniform ball
real_t relaxed_radius(index_t elements_per_cell)
{
	const simulation_parameters shape;
	const cell_type type = scale_to_elements(paper_cell, elements_per_cell, shape);

	// a twentieth of the largest stable timestep of an isolated pair
	const real_t pair_stiffness = 8.0 * type.intra_well_depth * shape.steepness * shape.steepness
								  / (type.equilibrium_distance * type.equilibrium_distance);
	const real_t timestep = 0.05 * type.drag_coefficient / pair_stiffness;

	environment env = make_environment(3, timestep);
	env.cell_types[0] = type;
	env.create_cell();

	// elements scattered uniformly in a ball of the requested radius
	std::mt19937 generator(1);
	std::uniform_real_distribution<real_t> coordinate(-1.0, 1.0);
	while (env.agents->size() < elements_per_cell)
	{
		const real_t x = coordinate(generator);
		const real_t y = coordinate(generator);
		const real_t z = coordinate(generator);
		if (x * x + y * y + z * z <= 1.0)
		{
			add_agent(env, 0, { paper_cell.radius * x, paper_cell.radius * y, paper_cell.radius * z });
		}
	}

	// four seconds, four times the relaxation time damping / stiffness of the cell
	openmp_solver solver;
	solver.solve(env, static_cast<index_t>(4.0 / timestep));

	const auto centre = centre_of_mass(env);
	real_t squared_gyration_radius = 0.0;
	for (index_t i = 0; i < env.agents->size(); ++i)
	{
		const auto position = env.agents->get_agent_at(i)->position();
		for (std::size_t d = 0; d < 3; ++d)
		{
			squared_gyration_radius += std::pow(position[d] - centre.at(d), 2) / static_cast<real_t>(elements_per_cell);
		}
	}
	return std::sqrt(5.0 / 3.0 * squared_gyration_radius);
}
} // namespace

TEST(ScaledCellTest, RelaxedCellHasTheRequestedRadiusForDifferentNumbersOfElements)
{
	const real_t coarse = relaxed_radius(128);
	const real_t fine = relaxed_radius(256);

	// With steepness 2 the elements compress each other slightly (as the paper notes), so the cell
	// comes out a little smaller than requested: about 9.3 for a requested radius of 10.
	EXPECT_NEAR(coarse, paper_cell.radius, 1.0);
	EXPECT_NEAR(fine, paper_cell.radius, 1.0);
	EXPECT_LT(coarse, paper_cell.radius);
	EXPECT_LT(fine, paper_cell.radius);

	// Doubling the number of elements changes the size of the cell by well under one percent of it,
	// although the element spacing shrinks by a fifth.
	EXPECT_NEAR(coarse, fine, 0.2);
}

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver::tests
