#include "micromechanics/agent.h"

#include <memory>
#include <tuple>

namespace physicore::mechanics::micromechanics {

agent::agent(index_t index, agent_data& data)
	: base_agent_interface(index),
	  physicore::base_agent_generic_storage<base_agent_data>(index, data.base_data),
	  data(data)
{}

agent::agent(index_t index, std::tuple<std::unique_ptr<base_agent_data>, std::unique_ptr<agent_data>>& datas)
	: agent(index, *std::get<std::unique_ptr<agent_data>>(datas))
{}

index_t& agent::cell_index() { return data.cell_indices[this->index]; }

std::span<real_t> agent::velocity()
{
	const index_t dims = data.base_data.dims;
	return { &data.velocities[this->index * dims], dims };
}

std::span<real_t> agent::previous_velocity()
{
	const index_t dims = data.base_data.dims;
	return { &data.previous_velocities[this->index * dims], dims };
}

} // namespace physicore::mechanics::micromechanics
