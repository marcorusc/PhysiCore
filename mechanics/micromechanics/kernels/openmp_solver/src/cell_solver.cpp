#include "cell_solver.h"

#include <algorithm>

#include "parallel_threshold.h"

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

void cell_solver::update(environment& e)
{
	const agent_data& agents = generic_agent_solver<agent>::retrieve_agent_data(*e.agents);
	cell_data& cells = generic_agent_solver<cell>::retrieve_agent_data(*e.cells);

	const index_t dims = cells.dims;
	const index_t cells_count = cells.cells_count;

	std::ranges::fill(cells.agent_counts, 0);
	for (index_t i = 0; i < agents.agents_count; ++i)
	{
		const index_t cell_index = agents.cell_indices[i];
		if (cell_index < cells_count)
		{
			++cells.agent_counts[cell_index];
		}
	}

	for (index_t c = 0; c < cells_count; ++c)
	{
		if (cells.agent_counts[c] != 0)
		{
			std::fill_n(&cells.positions[c * dims], dims, 0.0);
			std::fill_n(&cells.velocities[c * dims], dims, 0.0);
		}
	}

	// Summed by one thread in sub-agent order, so the result does not depend on the number of threads.
	for (index_t i = 0; i < agents.agents_count; ++i)
	{
		const index_t cell_index = agents.cell_indices[i];
		if (cell_index >= cells_count)
			continue;

		for (index_t d = 0; d < dims; ++d)
		{
			cells.positions[cell_index * dims + d] += agents.base_data.positions[i * dims + d];
			cells.velocities[cell_index * dims + d] += agents.velocities[i * dims + d];
		}
	}

#pragma omp parallel for if (cells_count >= parallel_threshold)
	for (index_t c = 0; c < cells_count; ++c)
	{
		if (cells.agent_counts[c] == 0)
			continue;

		const real_t inverse_count = 1.0 / static_cast<real_t>(cells.agent_counts[c]);
		for (index_t d = 0; d < dims; ++d)
		{
			cells.positions[c * dims + d] *= inverse_count;
			cells.velocities[c * dims + d] *= inverse_count;
		}
	}
}

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
