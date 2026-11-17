#include "micromechanics/agent_data.h"

#include <cassert>

#include "micromechanics/agent_interface.h"

namespace physicore::mechanics::micromechanics {

agent_data::agent_data(physicore::base_agent_data& base_data) : base_data(base_data) {}

void agent_data::add()
{
	agents_count = base_data.agents_count;
	const index_t dims = base_data.dims;

	cell_indices.resize(agents_count, no_cell);
	velocities.resize(agents_count * dims, 0.0);
}

void agent_data::remove_at(index_t position)
{
	assert(position < agents_count);
	if (position >= agents_count)
		return;

	const index_t dims = base_data.dims;
	const index_t last = agents_count - 1;

	if (position != last)
	{
		physicore::base_agent_data::move_scalar(&cell_indices[position], &cell_indices[last]);
		physicore::base_agent_data::move_vector(&velocities[position * dims], &velocities[last * dims], dims);
	}

	agents_count = last;

	cell_indices.resize(agents_count);
	velocities.resize(agents_count * dims);
}

} // namespace physicore::mechanics::micromechanics
