#include "position_solver.h"

#include <cassert>
#include <cstdint>

#include "parallel_threshold.h"

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

void position_solver::update(environment& e, std::span<const real_t> forces)
{
	agent_data& data = generic_agent_solver<agent>::retrieve_agent_data(*e.agents);
	const cell_data& cells = generic_agent_solver<cell>::retrieve_agent_data(*e.cells);

	const index_t agents_count = data.agents_count;
	const index_t dims = data.base_data.dims;

	assert(forces.size() == agents_count * dims);

	const real_t inverse_drag = 1.0 / e.params.drag_coefficient;
	const real_t timestep = e.mechanics_timestep;

	real_t* positions = data.base_data.positions.data();
	real_t* velocities = data.velocities.data();
	real_t* previous_velocities = data.previous_velocities.data();
	const index_t* cell_indices = data.cell_indices.data();
	const std::uint8_t* is_movable = cells.is_movable.data();
	const real_t* force = forces.data();

#pragma omp parallel for if (agents_count >= parallel_threshold)
	for (index_t i = 0; i < agents_count; ++i)
	{
		const index_t cell_index = cell_indices[i];
		const bool movable = cell_index == no_cell || is_movable[cell_index] != 0;

		for (index_t d = 0; d < dims; ++d)
		{
			const index_t k = i * dims + d;

			previous_velocities[k] = velocities[k];
			velocities[k] = movable ? force[k] * inverse_drag : 0.0;
			positions[k] += velocities[k] * timestep;
		}
	}
}

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
