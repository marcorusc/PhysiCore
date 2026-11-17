#include "pressure_solver.h"

#include <algorithm>
#include <numbers>

#include <micromechanics/pair_potential.h>

#include "parallel_threshold.h"

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

namespace {
// Everything needed for the contribution of one pair of sub-agents to the virial
struct outside_forces
{
	const real_t* positions;
	const index_t* cell_indices;
	const real_t* cell_positions;
	index_t dims;
	pair_potential inter;

	// (x_i - X_c) . F_ij for sub-agent i of cell c and a sub-agent j that is not in that cell
	real_t virial(index_t i, index_t j) const
	{
		const index_t c = cell_indices[i];

		real_t squared_distance = 0.0;
		real_t projection = 0.0;
		for (index_t d = 0; d < dims; ++d)
		{
			const real_t difference = positions[i * dims + d] - positions[j * dims + d];
			squared_distance += difference * difference;
			projection += (positions[i * dims + d] - cell_positions[c * dims + d]) * difference;
		}

		// the force on i is force_coefficient * (x_i - x_j)
		return force_coefficient(inter, squared_distance) * projection;
	}
};

// Volume that one sub-agent takes in a close packing at the equilibrium distance
real_t volume_per_element(index_t dims, real_t equilibrium_distance)
{
	switch (dims)
	{
		case 1:
			return equilibrium_distance;
		case 2:
			return 0.5 * std::numbers::sqrt3 * equilibrium_distance * equilibrium_distance;
		default:
			return equilibrium_distance * equilibrium_distance * equilibrium_distance / std::numbers::sqrt2;
	}
}
} // namespace

void pressure_solver::update(environment& e)
{
	const agent_data& agents = generic_agent_solver<agent>::retrieve_agent_data(*e.agents);
	const cell_data& cells = generic_agent_solver<cell>::retrieve_agent_data(*e.cells);

	const index_t agents_count = agents.agents_count;
	const index_t dims = agents.base_data.dims;

	const outside_forces pairs { .positions = agents.base_data.positions.data(),
								 .cell_indices = agents.cell_indices.data(),
								 .cell_positions = cells.positions.data(),
								 .dims = dims,
								 .inter = e.params.inter_potential() };

	grid_.build(pairs.positions, agents_count, dims, pairs.inter.cutoff_distance);

	virials_.resize(agents_count);
	real_t* virials = virials_.data();
	const neighbor_grid& grid = grid_;
	const index_t cells_count = cells.cells_count;

#pragma omp parallel for if (agents_count >= parallel_threshold)
	for (index_t i = 0; i < agents_count; ++i)
	{
		real_t virial = 0.0;

		if (pairs.cell_indices[i] < cells_count)
		{
			grid.visit_candidates(pairs.positions + i * dims, [&pairs, &virial, i](index_t j) {
				if (pairs.cell_indices[i] != pairs.cell_indices[j])
				{
					virial += pairs.virial(i, j);
				}
			});
		}

		virials[i] = virial;
	}

	accumulate(e);
}

void pressure_solver::update_brute_force(environment& e)
{
	const agent_data& agents = generic_agent_solver<agent>::retrieve_agent_data(*e.agents);
	const cell_data& cells = generic_agent_solver<cell>::retrieve_agent_data(*e.cells);

	const index_t agents_count = agents.agents_count;

	const outside_forces pairs { .positions = agents.base_data.positions.data(),
								 .cell_indices = agents.cell_indices.data(),
								 .cell_positions = cells.positions.data(),
								 .dims = agents.base_data.dims,
								 .inter = e.params.inter_potential() };

	virials_.assign(agents_count, 0.0);

	for (index_t i = 0; i < agents_count; ++i)
	{
		if (pairs.cell_indices[i] >= cells.cells_count)
		{
			continue;
		}
		for (index_t j = 0; j < agents_count; ++j)
		{
			if (pairs.cell_indices[i] != pairs.cell_indices[j])
			{
				virials_[i] += pairs.virial(i, j);
			}
		}
	}

	accumulate(e);
}

void pressure_solver::accumulate(environment& e)
{
	const agent_data& agents = generic_agent_solver<agent>::retrieve_agent_data(*e.agents);
	cell_data& cells = generic_agent_solver<cell>::retrieve_agent_data(*e.cells);

	const index_t dims = agents.base_data.dims;
	const index_t cells_count = cells.cells_count;

	std::ranges::fill(cells.pressures, 0.0);

	// Summed by one thread in sub-agent order, so the result does not depend on the number of threads.
	// The number of sub-agents of every cell is counted here and not taken from the cell, so the
	// pressure is right even if the stored count is out of date.
	std::vector<index_t> counts(cells_count, 0);
	for (index_t i = 0; i < agents.agents_count; ++i)
	{
		const index_t c = agents.cell_indices[i];
		if (c < cells_count)
		{
			cells.pressures[c] += virials_[i];
			++counts[c];
		}
	}

	const real_t element_volume = volume_per_element(dims, e.params.equilibrium_distance);

	for (index_t c = 0; c < cells_count; ++c)
	{
		if (counts[c] == 0)
		{
			continue;
		}
		const real_t volume = static_cast<real_t>(counts[c]) * element_volume;
		cells.pressures[c] = -cells.pressures[c] / (static_cast<real_t>(dims) * volume);
	}
}

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
