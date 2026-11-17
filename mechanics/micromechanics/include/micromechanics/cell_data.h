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

	std::vector<std::uint8_t> is_movable;
	std::vector<std::uint8_t> is_motile;

	std::vector<real_t> migration_speeds;
	std::vector<real_t> migration_biases;
	std::vector<real_t> migration_bias_directions;
	std::vector<real_t> persistence_times;
	std::vector<real_t> motility_directions;

	std::vector<real_t> positions;
	std::vector<real_t> velocities;
	std::vector<index_t> agent_counts;
	std::vector<real_t> pressures;

	std::vector<std::vector<index_t>> neighbor_cells;
};

} // namespace physicore::mechanics::micromechanics
