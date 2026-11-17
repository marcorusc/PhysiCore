#include "micromechanics/environment.h"

#include <algorithm>
#include <memory>
#include <stdexcept>

#include <common/base_agent_data.h>

#include "micromechanics/cell_placement.h"

namespace physicore::mechanics::micromechanics {

environment::environment(const cartesian_mesh& mesh, real_t timestep) : mechanics_timestep(timestep), mesh(mesh)
{
	cells = std::make_shared<cell_container>(std::make_unique<cell_data>(mesh.dims));

	auto base_data = std::make_unique<physicore::base_agent_data>(mesh.dims);
	auto data = std::make_unique<agent_data>(*base_data);
	agents = std::make_shared<agent_container>(std::move(base_data), std::move(data));
}

void environment::run_single_timestep()
{
	if (!solver)
	{
		throw std::runtime_error("Micromechanics solver is not set");
	}

	solver->solve(*this, 1);
}

void environment::serialize_state(real_t current_time)
{
	if (serializer)
	{
		serializer->serialize(*this, current_time);
	}
}

cell_interface* environment::create_cell()
{
	auto* new_cell = cells->create();
	new_cell->agent_radius() = 0.5 * params.equilibrium_distance;
	return new_cell;
}

cell_interface* environment::create_cell(std::span<const real_t> centre, index_t elements_count)
{
	const index_t dims = mesh.dims;
	if (centre.size() != dims)
	{
		throw std::invalid_argument("The centre of a cell needs one value per dimension");
	}

	const std::vector<real_t> offsets = packed_ball_positions(dims, elements_count, params.equilibrium_distance);

	auto* new_cell = create_cell();
	const index_t cell_index = cells->size() - 1;

	// the cell sits at the centre of mass of its sub-agents
	std::ranges::copy(centre, new_cell->position().begin());
	for (index_t i = 0; i < elements_count; ++i)
	{
		auto position = create_agent(cell_index)->position();
		for (index_t d = 0; d < dims; ++d)
		{
			position[d] = centre[d] + offsets[i * dims + d];
			new_cell->position()[d] += offsets[i * dims + d] / static_cast<real_t>(elements_count);
		}
	}

	return new_cell;
}

agent_interface* environment::create_agent(index_t cell_index)
{
	if (cell_index >= cells->size())
	{
		throw std::out_of_range("Invalid cell index");
	}

	auto* new_agent = agents->create();
	new_agent->cell_index() = cell_index;

	++generic_agent_solver<cell>::retrieve_agent_data(*cells).agent_counts[cell_index];

	return new_agent;
}

void environment::remove_agent(index_t agent_index)
{
	if (agent_index >= agents->size())
	{
		throw std::out_of_range("Invalid agent index");
	}

	const index_t cell_index = generic_agent_solver<agent>::retrieve_agent_data(*agents).cell_indices[agent_index];
	if (cell_index != no_cell)
	{
		--generic_agent_solver<cell>::retrieve_agent_data(*cells).agent_counts[cell_index];
	}

	agents->remove_at(agent_index);
}

void environment::remove_cell(index_t cell_index)
{
	if (cell_index >= cells->size())
	{
		throw std::out_of_range("Invalid cell index");
	}

	auto& cell_indices = generic_agent_solver<agent>::retrieve_agent_data(*agents).cell_indices;

	// Going backwards, the sub-agent moved into a freed slot has already been checked.
	for (index_t i = agents->size(); i-- > 0;)
	{
		if (cell_indices[i] == cell_index)
		{
			agents->remove_at(i);
		}
	}

	// The last cell is about to take over the freed index.
	const index_t last_cell_index = cells->size() - 1;
	if (cell_index != last_cell_index)
	{
		std::ranges::replace(cell_indices, last_cell_index, cell_index);
	}

	cells->remove_at(cell_index);
}

} // namespace physicore::mechanics::micromechanics
