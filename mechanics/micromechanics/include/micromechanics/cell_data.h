#pragma once

#include <cstdint>
#include <vector>

#include <common/types.h>

namespace physicore::mechanics::micromechanics {

struct cell_data
{
	explicit cell_data(index_t dims);

	void add();
	void remove_at(index_t position);

	index_t cells_count = 0;
	index_t dims;

	std::vector<index_t> type_indices;
	std::vector<std::uint8_t> is_movable;

	std::vector<real_t> positions;
	std::vector<real_t> velocities;
	std::vector<index_t> agent_counts;

	std::vector<std::vector<index_t>> neighbor_cells;
};

} // namespace physicore::mechanics::micromechanics
