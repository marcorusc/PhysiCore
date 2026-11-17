#include "micromechanics/cell_data.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <utility>

namespace physicore::mechanics::micromechanics {

cell_data::cell_data(index_t dims) : dims(dims) { assert(dims > 0); }

void cell_data::add()
{
	++cells_count;
	const std::size_t count = static_cast<std::size_t>(cells_count);
	const std::size_t vector_size = count * static_cast<std::size_t>(dims);

	agent_radius.resize(count, 1.0);
	is_movable.resize(count, 1);
	is_motile.resize(count, 0);

	migration_speeds.resize(count, 0.0);
	migration_biases.resize(count, 0.0);
	migration_bias_directions.resize(vector_size, 0.0);
	persistence_times.resize(count, 0.0);
	motility_directions.resize(vector_size, 0.0);

	positions.resize(vector_size, 0.0);
	velocities.resize(vector_size, 0.0);
	agent_counts.resize(count, 0);
	pressures.resize(count, 0.0);

	neighbor_cells.resize(count);
}

void cell_data::remove_at(index_t position)
{
	assert(position < cells_count);
	if (position >= cells_count)
		return;

	const index_t last = cells_count - 1;
	if (position != last)
	{
		agent_radius[position] = agent_radius[last];
		is_movable[position] = is_movable[last];
		is_motile[position] = is_motile[last];

		migration_speeds[position] = migration_speeds[last];
		migration_biases[position] = migration_biases[last];
		persistence_times[position] = persistence_times[last];
		agent_counts[position] = agent_counts[last];
		pressures[position] = pressures[last];

		const auto width = static_cast<std::size_t>(dims);
		const auto destination = static_cast<std::size_t>(position) * width;
		const auto source = static_cast<std::size_t>(last) * width;
		std::copy_n(migration_bias_directions.begin() + source, width, migration_bias_directions.begin() + destination);
		std::copy_n(motility_directions.begin() + source, width, motility_directions.begin() + destination);
		std::copy_n(positions.begin() + source, width, positions.begin() + destination);
		std::copy_n(velocities.begin() + source, width, velocities.begin() + destination);

		neighbor_cells[position] = std::move(neighbor_cells[last]);
	}

	cells_count = last;
	const std::size_t count = static_cast<std::size_t>(cells_count);
	const std::size_t vector_size = count * static_cast<std::size_t>(dims);

	agent_radius.resize(count);
	is_movable.resize(count);
	is_motile.resize(count);

	migration_speeds.resize(count);
	migration_biases.resize(count);
	migration_bias_directions.resize(vector_size);
	persistence_times.resize(count);
	motility_directions.resize(vector_size);

	positions.resize(vector_size);
	velocities.resize(vector_size);
	agent_counts.resize(count);
	pressures.resize(count);

	neighbor_cells.resize(count);
}

} // namespace physicore::mechanics::micromechanics
