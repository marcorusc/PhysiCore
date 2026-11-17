#include "micromechanics/cell.h"

#include <cassert>
#include <memory>
#include <tuple>

namespace physicore::mechanics::micromechanics {

cell::cell(index_t index, cell_data& data) : base_agent_interface(index), data(data)
{
	assert(index < data.cells_count);
}

cell::cell(index_t index, std::tuple<std::unique_ptr<cell_data>>& datas)
	: cell(index, *std::get<std::unique_ptr<cell_data>>(datas))
{}

std::span<real_t> cell::position()
{
	return { &data.positions[this->index * data.dims], static_cast<std::size_t>(data.dims) };
}

real_t& cell::agent_radius() { return data.agent_radius[this->index]; }

std::uint8_t& cell::is_movable() { return data.is_movable[this->index]; }

std::uint8_t& cell::is_motile() { return data.is_motile[this->index]; }

real_t& cell::migration_speed() { return data.migration_speeds[this->index]; }

real_t& cell::migration_bias() { return data.migration_biases[this->index]; }

std::span<real_t> cell::migration_bias_direction()
{
	return { &data.migration_bias_directions[this->index * data.dims], static_cast<std::size_t>(data.dims) };
}

real_t& cell::persistence_time() { return data.persistence_times[this->index]; }

std::span<real_t> cell::motility_direction()
{
	return { &data.motility_directions[this->index * data.dims], static_cast<std::size_t>(data.dims) };
}

std::span<real_t> cell::velocity()
{
	return { &data.velocities[this->index * data.dims], static_cast<std::size_t>(data.dims) };
}

index_t cell::agent_count() const { return data.agent_counts[this->index]; }

real_t cell::pressure() const { return data.pressures[this->index]; }

std::span<index_t> cell::neighbors() { return { data.neighbor_cells[this->index] }; }

} // namespace physicore::mechanics::micromechanics
