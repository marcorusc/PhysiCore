#include "micromechanics/cell_data.h"

#include <cassert>
#include <utility>

#include <common/base_agent_data.h>

namespace physicore::mechanics::micromechanics {

cell_data::cell_data(index_t dims) : dims(dims) { assert(dims > 0); }

void cell_data::add()
{
	++cells_count;

	type_indices.resize(cells_count, 0);
	is_movable.resize(cells_count, 1);

	positions.resize(cells_count * dims, 0.0);
	velocities.resize(cells_count * dims, 0.0);
	agent_counts.resize(cells_count, 0);

	neighbor_cells.resize(cells_count);
}

void cell_data::remove_at(index_t position)
{
	assert(position < cells_count);
	if (position >= cells_count)
		return;

	using storage = physicore::base_agent_data;

	const index_t last = cells_count - 1;
	if (position != last)
	{
		storage::move_scalar(&type_indices[position], &type_indices[last]);
		storage::move_scalar(&is_movable[position], &is_movable[last]);

		storage::move_vector(&positions[position * dims], &positions[last * dims], dims);
		storage::move_vector(&velocities[position * dims], &velocities[last * dims], dims);
		storage::move_scalar(&agent_counts[position], &agent_counts[last]);

		neighbor_cells[position] = std::move(neighbor_cells[last]);
	}

	cells_count = last;

	type_indices.resize(cells_count);
	is_movable.resize(cells_count);

	positions.resize(cells_count * dims);
	velocities.resize(cells_count * dims);
	agent_counts.resize(cells_count);

	neighbor_cells.resize(cells_count);
}

} // namespace physicore::mechanics::micromechanics
