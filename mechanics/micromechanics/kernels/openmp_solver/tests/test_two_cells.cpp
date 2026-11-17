#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

#include <gtest/gtest.h>

#include "openmp_solver.h"
#include "test_helpers.h"

// Two cells of 27 sub-agents each, placed side by side along x with a gap that is inside the interaction
// range. What happens next is decided by the ratio of the inter- to the intra-cellular well depth.

namespace physicore::mechanics::micromechanics::kernels::openmp_solver::tests {

namespace {
constexpr real_t initial_gap = 1.3;

void add_two_cells(environment& env)
{
	add_jittered_lattice_cell(env, 1.0, 0.15, 0.0, 42);
	add_jittered_lattice_cell(env, 1.0, 0.15, 2.0 + initial_gap, 43);
}

struct two_cell_state
{
	// distance between the centres of mass of the two cells
	real_t centre_distance = 0.0;
	// smallest x of cell 1 minus largest x of cell 0: positive while no element has passed an element
	// of the other cell along the axis joining the cells
	real_t gap_along_axis = 0.0;
	// number of pairs from different cells that are closer than 1.3 equilibrium distances
	index_t contacts = 0;
	// largest distance of an element of cell 0 from the centre of cell 0
	real_t radius_of_first_cell = 0.0;
	real_t largest_speed = 0.0;
};

two_cell_state measure(environment& env)
{
	std::array<std::array<real_t, 3>, 2> centres {};
	std::array<index_t, 2> counts {};
	real_t largest_x_of_first = std::numeric_limits<real_t>::lowest();
	real_t smallest_x_of_second = std::numeric_limits<real_t>::max();

	two_cell_state state;
	for (index_t i = 0; i < env.agents->size(); ++i)
	{
		auto* element = env.agents->get_agent_at(i);
		const index_t cell_index = element->cell_index();
		const auto position = element->position();

		++counts.at(cell_index);
		for (std::size_t d = 0; d < 3; ++d)
		{
			centres.at(cell_index).at(d) += position[d];
			state.largest_speed = std::max(state.largest_speed, std::abs(element->velocity()[d]));
		}
		if (cell_index == 0)
			largest_x_of_first = std::max(largest_x_of_first, position[0]);
		else
			smallest_x_of_second = std::min(smallest_x_of_second, position[0]);
	}
	state.gap_along_axis = smallest_x_of_second - largest_x_of_first;

	real_t squared_centre_distance = 0.0;
	for (std::size_t d = 0; d < 3; ++d)
	{
		centres[0].at(d) /= static_cast<real_t>(counts[0]);
		centres[1].at(d) /= static_cast<real_t>(counts[1]);
		squared_centre_distance += std::pow(centres[0].at(d) - centres[1].at(d), 2);
	}
	state.centre_distance = std::sqrt(squared_centre_distance);

	for (index_t i = 0; i < env.agents->size(); ++i)
	{
		auto* element = env.agents->get_agent_at(i);
		if (element->cell_index() == 0)
		{
			real_t squared_radius = 0.0;
			for (std::size_t d = 0; d < 3; ++d)
			{
				squared_radius += std::pow(element->position()[d] - centres[0].at(d), 2);
			}
			state.radius_of_first_cell = std::max(state.radius_of_first_cell, std::sqrt(squared_radius));
		}

		for (index_t j = i + 1; j < env.agents->size(); ++j)
		{
			if (element->cell_index() != env.agents->get_agent_at(j)->cell_index()
				&& squared_distance(env, i, j) < 1.3 * 1.3)
			{
				++state.contacts;
			}
		}
	}

	return state;
}
} // namespace

TEST(TwoCellsTest, DefaultParametersMakeCellsAdhereWithoutMixing)
{
	environment env = make_environment(3, 0.002);
	ASSERT_DOUBLE_EQ(env.params.inter_well_depth / env.params.intra_well_depth, 0.5);
	add_two_cells(env);
	const two_cell_state before = measure(env);

	openmp_solver solver;
	solver.solve(env, 5000);

	const two_cell_state after = measure(env);

	// the cells pull each other closer and build up contacts
	EXPECT_LT(after.centre_distance, before.centre_distance - 0.5);
	EXPECT_GT(after.contacts, 3 * before.contacts);

	// but stay two cells: no element has passed an element of the other cell,
	EXPECT_GT(after.gap_along_axis, 0.1);
	// and each cell is still as compact as an isolated one (27 elements fill a ball of radius about 1.7)
	EXPECT_LT(after.radius_of_first_cell, 2.2);

	// and the pair has come to rest
	EXPECT_LT(after.largest_speed, 1e-3);
}

TEST(TwoCellsTest, WithoutInterCellularWellCellsIgnoreEachOther)
{
	// The well depth scales repulsion as well as adhesion, so a zero inter-cellular depth removes the
	// interaction between cells altogether: each cell evolves exactly as if it were alone.
	environment both = make_environment(3, 0.002);
	both.params.inter_well_depth = 0.0;
	add_two_cells(both);

	environment alone = make_environment(3, 0.002);
	add_jittered_lattice_cell(alone, 1.0, 0.15, 0.0, 42);

	const real_t initial_centre_distance = measure(both).centre_distance;

	openmp_solver solver_for_both;
	openmp_solver solver_for_one;
	solver_for_both.solve(both, 500);
	solver_for_one.solve(alone, 500);

	EXPECT_NEAR(measure(both).centre_distance, initial_centre_distance, 1e-10);

	const auto positions_of_both = positions_of(both);
	const auto positions_alone = positions_of(alone);
	for (std::size_t k = 0; k < positions_alone.size(); ++k)
	{
		EXPECT_NEAR(positions_of_both[k], positions_alone[k], 1e-10);
	}
}

TEST(TwoCellsTest, EqualWellDepthsMakeCellsIndistinguishableAndTheyMerge)
{
	// With equal depths every pair interacts in the same way, whatever cell the elements belong to...
	environment two_cells = make_environment(3, 0.002);
	two_cells.params.inter_well_depth = two_cells.params.intra_well_depth;
	add_two_cells(two_cells);

	// ... so the dynamics is that of a single cell holding all the elements.
	environment one_cell = make_environment(3, 0.002);
	add_two_cells(one_cell);
	for (index_t i = 0; i < one_cell.agents->size(); ++i)
	{
		one_cell.agents->get_agent_at(i)->cell_index() = 0;
	}

	openmp_solver solver_for_two;
	openmp_solver solver_for_one;
	solver_for_two.solve(two_cells, 5000);
	solver_for_one.solve(one_cell, 5000);

	EXPECT_EQ(positions_of(two_cells), positions_of(one_cell));

	// Nothing keeps the cells apart: elements of one cell move in between elements of the other.
	EXPECT_LT(measure(two_cells).gap_along_axis, -0.1);
}

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver::tests
