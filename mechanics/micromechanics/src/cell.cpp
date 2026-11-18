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

index_t& cell::type_index() { return data.type_indices[this->index]; }

std::uint8_t& cell::is_movable() { return data.is_movable[this->index]; }

std::span<real_t> cell::velocity()
{
	return { &data.velocities[this->index * data.dims], static_cast<std::size_t>(data.dims) };
}

index_t cell::agent_count() const { return data.agent_counts[this->index]; }

std::span<index_t> cell::neighbors() { return { data.neighbor_cells[this->index] }; }

} // namespace physicore::mechanics::micromechanics
