#include "position_solver.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>

#include <micromechanics/counter_random.h>

#include "parallel_threshold.h"

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

void position_solver::update(environment& e, std::span<const real_t> forces)
{
	agent_data& data = generic_agent_solver<agent>::retrieve_agent_data(*e.agents);
	const cell_data& cells = generic_agent_solver<cell>::retrieve_agent_data(*e.cells);

	const index_t agents_count = data.agents_count;
	const index_t dims = data.base_data.dims;

	assert(forces.size() == agents_count * dims);

	const real_t inverse_drag = 1.0 / e.params.drag_coefficient;
	const real_t timestep = e.mechanics_timestep;

	const bool noisy = e.params.diffusion_coefficient > 0.0;
	const real_t noise_amplitude = noisy ? std::sqrt(2.0 * e.params.diffusion_coefficient * timestep) : 0.0;
	const std::uint64_t seed = e.params.random_seed;
	const index_t step = step_;

	real_t* positions = data.base_data.positions.data();
	real_t* velocities = data.velocities.data();
	const index_t* cell_indices = data.cell_indices.data();
	const std::uint8_t* is_movable = cells.is_movable.data();
	const std::uint8_t* is_motile = cells.is_motile.data();
	const real_t* migration_speeds = cells.migration_speeds.data();
	const real_t* motility_directions = cells.motility_directions.data();
	const real_t* force = forces.data();

	const bool walls = e.params.domain_walls;
	const cartesian_mesh& mesh = e.mesh;
	index_t stops = 0;

#pragma omp parallel for reduction(+ : stops) if (agents_count >= parallel_threshold)
	for (index_t i = 0; i < agents_count; ++i)
	{
		const index_t cell_index = cell_indices[i];
		const bool movable = cell_index == no_cell || is_movable[cell_index] != 0;
		const bool motile = cell_index != no_cell && is_motile[cell_index] != 0;

		for (index_t d = 0; d < dims; ++d)
		{
			const index_t k = i * dims + d;

			const real_t crawling =
				motile ? migration_speeds[cell_index] * motility_directions[cell_index * dims + d] : 0.0;

			velocities[k] = movable ? force[k] * inverse_drag + crawling : 0.0;
			positions[k] += velocities[k] * timestep;

			if (noisy && movable)
			{
				positions[k] += noise_amplitude * standard_normal(seed, step, i, d);
			}

			if (walls && movable)
			{
				const real_t inside = std::clamp(positions[k], static_cast<real_t>(mesh.bounding_box_mins.at(d)),
												 static_cast<real_t>(mesh.bounding_box_maxs.at(d)));
				if (inside != positions[k])
				{
					positions[k] = inside;
					++stops;
				}
			}
		}
	}

	wall_stops_ += stops;
	++step_;
}

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
