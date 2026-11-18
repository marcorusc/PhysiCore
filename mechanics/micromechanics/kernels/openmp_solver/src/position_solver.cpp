#include "position_solver.h"

#include <cassert>
#include <cstdint>
#include <vector>

#include "element_types.h"
#include "parallel_threshold.h"

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

void position_solver::update(environment& e, std::span<const real_t> forces)
{
	agent_data& data = generic_agent_solver<agent>::retrieve_agent_data(*e.agents);
	const cell_data& cells = generic_agent_solver<cell>::retrieve_agent_data(*e.cells);

	const index_t agents_count = data.agents_count;
	const index_t dims = data.base_data.dims;

	assert(forces.size() == agents_count * dims);

	const real_t timestep = e.mechanics_timestep;

	// per cell type: one over the drag of an element
	std::vector<real_t> inverse_drags;
	inverse_drags.reserve(e.cell_types.size());
	for (const cell_type& type : e.cell_types)
	{
		inverse_drags.push_back(1.0 / type.drag_coefficient);
	}
	const real_t* inverse_drag_of_type = inverse_drags.data();
	const index_t* type_indices = cells.type_indices.data();
	const index_t cells_count = cells.cells_count;

	real_t* positions = data.base_data.positions.data();
	real_t* velocities = data.velocities.data();
	const index_t* cell_indices = data.cell_indices.data();
	const std::uint8_t* is_movable = cells.is_movable.data();
	const real_t* force = forces.data();

#pragma omp parallel for if (agents_count >= parallel_threshold)
	for (index_t i = 0; i < agents_count; ++i)
	{
		const index_t cell_index = cell_indices[i];
		const bool movable = cell_index == no_cell || is_movable[cell_index] != 0;
		const index_t type = type_of_element(cell_index, type_indices, cells_count);
		const real_t inverse_drag = inverse_drag_of_type[type];

		for (index_t d = 0; d < dims; ++d)
		{
			const index_t k = i * dims + d;

			velocities[k] = movable ? force[k] * inverse_drag : 0.0;
			positions[k] += velocities[k] * timestep;
		}
	}
}

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
