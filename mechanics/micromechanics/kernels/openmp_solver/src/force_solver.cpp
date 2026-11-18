#include "force_solver.h"

#include <array>

#include <micromechanics/cell_type.h>
#include <micromechanics/pair_potential.h>

#include "element_types.h"
#include "parallel_threshold.h"

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

namespace {
// Everything needed to evaluate the force between two sub-agents
struct pair_forces
{
	const real_t* positions;
	const index_t* cell_indices;
	const index_t* type_indices;
	index_t cells_count;
	index_t dims;
	const interaction_table* table;

	// Adds the force on sub-agent i due to sub-agent j to net_force
	void add(index_t i, index_t j, std::array<real_t, 3>& net_force) const
	{
		std::array<real_t, 3> difference {};
		real_t squared_distance = 0.0;
		for (index_t d = 0; d < dims; ++d)
		{
			difference[d] = positions[i * dims + d] - positions[j * dims + d];
			squared_distance += difference[d] * difference[d];
		}

		const index_t type_i = type_of_element(cell_indices[i], type_indices, cells_count);
		const bool same_cell = cell_indices[i] == cell_indices[j] && cell_indices[i] != no_cell;

		const pair_interaction& interaction =
			same_cell ? table->intra(type_i)
					  : table->inter(type_i, type_of_element(cell_indices[j], type_indices, cells_count));
		const real_t coefficient = force_coefficient(interaction, squared_distance);

		for (index_t d = 0; d < dims; ++d)
		{
			net_force[d] += coefficient * difference[d];
		}
	}
};

pair_forces make_pair_forces(const agent_data& data, const cell_data& cells, const interaction_table& table)
{
	return { .positions = data.base_data.positions.data(),
			 .cell_indices = data.cell_indices.data(),
			 .type_indices = cells.type_indices.data(),
			 .cells_count = cells.cells_count,
			 .dims = data.base_data.dims,
			 .table = &table };
}
} // namespace

void force_solver::compute(environment& e)
{
	const agent_data& data = generic_agent_solver<agent>::retrieve_agent_data(*e.agents);
	const cell_data& cells = generic_agent_solver<cell>::retrieve_agent_data(*e.cells);
	const interaction_table table(e.cell_types, e.params);
	const pair_forces pairs = make_pair_forces(data, cells, table);

	const index_t agents_count = data.agents_count;
	const index_t dims = pairs.dims;

	// bins as wide as the longest reach of any pair of cell types
	grid_.build(pairs.positions, agents_count, dims, table.largest_cutoff());

	forces_.resize(agents_count * dims);
	real_t* forces = forces_.data();
	const neighbor_grid& grid = grid_;

#pragma omp parallel for if (agents_count >= parallel_threshold)
	for (index_t i = 0; i < agents_count; ++i)
	{
		std::array<real_t, 3> net_force {};

		grid.visit_candidates(pairs.positions + i * dims, [&pairs, &net_force, i](index_t j) {
			if (i != j)
			{
				pairs.add(i, j, net_force);
			}
		});

		for (index_t d = 0; d < dims; ++d)
		{
			forces[i * dims + d] = net_force[d];
		}
	}
}

void force_solver::compute_brute_force(environment& e)
{
	const agent_data& data = generic_agent_solver<agent>::retrieve_agent_data(*e.agents);
	const cell_data& cells = generic_agent_solver<cell>::retrieve_agent_data(*e.cells);
	const interaction_table table(e.cell_types, e.params);
	const pair_forces pairs = make_pair_forces(data, cells, table);

	const index_t agents_count = data.agents_count;
	const index_t dims = pairs.dims;

	forces_.resize(agents_count * dims);
	real_t* forces = forces_.data();

#pragma omp parallel for if (agents_count >= parallel_threshold)
	for (index_t i = 0; i < agents_count; ++i)
	{
		std::array<real_t, 3> net_force {};

		for (index_t j = 0; j < agents_count; ++j)
		{
			if (i != j)
			{
				pairs.add(i, j, net_force);
			}
		}

		for (index_t d = 0; d < dims; ++d)
		{
			forces[i * dims + d] = net_force[d];
		}
	}
}

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
