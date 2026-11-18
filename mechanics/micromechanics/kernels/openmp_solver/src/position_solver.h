#pragma once

#include <span>

#include <common/generic_agent_solver.h>
#include <common/types.h>
#include <micromechanics/environment.h>

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

/**
 * @brief Position update of the overdamped Langevin equation (Euler-Maruyama scheme).
 *
 *     velocity = force / drag_coefficient + migration_speed * motility_direction
 *     position += velocity * timestep + sqrt(2 * thermal_energy / drag_coefficient * timestep) * xi
 *
 * with drag_coefficient that of the sub-agent's cell type, and
 * xi an independent standard normal number for every sub-agent, axis and step
 * (Sandersius and Newman 2008, Eq. 1 and 2). Without noise this is forward Euler. velocity is the
 * deterministic drift only; the noise has no finite velocity.
 *
 * The second term of the velocity is the motility. It is there only for sub-agents of a motile cell
 * and is the same for all of them, so on its own it moves the cell without changing its shape. It is
 * equivalent to a traction force drag_coefficient * migration_speed * motility_direction on every
 * sub-agent. The direction is chosen by motility_solver.
 *
 * Sub-agents of a cell that is not movable keep their position and get zero velocity.
 *
 * If the domain has walls, a sub-agent that ends the step outside the domain is put back on the face
 * it crossed. The wall force normally keeps sub-agents inside (see domain_wall.h); this is the guard
 * for a push or a noise kick stronger than the wall. wall_stops() counts how often it was needed.
 */
class position_solver : private generic_agent_solver<agent>, private generic_agent_solver<cell>
{
	// number of updates done so far; together with the seed it identifies the noise of a step
	index_t step_ = 0;

	index_t wall_stops_ = 0;

public:
	// Number of times a sub-agent was put back on a wall, summed over all updates and axes
	index_t wall_stops() const { return wall_stops_; }

	// forces holds agents_count * dims values, laid out like the positions
	void update(environment& e, std::span<const real_t> forces);
};

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
