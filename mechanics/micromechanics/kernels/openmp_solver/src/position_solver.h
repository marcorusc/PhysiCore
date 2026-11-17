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
 *     position += velocity * timestep + sqrt(2 * diffusion_coefficient * timestep) * xi
 *
 * with xi an independent standard normal number for every sub-agent, axis and step
 * (Sandersius and Newman 2008, Eq. 1 and 2). Without noise this is forward Euler. velocity is the
 * deterministic drift only; the noise has no finite velocity.
 *
 * The second term of the velocity is the motility. It is there only for sub-agents of a motile cell
 * and is the same for all of them, so on its own it moves the cell without changing its shape. It is
 * equivalent to a traction force drag_coefficient * migration_speed * motility_direction on every
 * sub-agent. The direction is chosen by motility_solver.
 *
 * Sub-agents of a cell that is not movable keep their position and get zero velocity.
 */
class position_solver : private generic_agent_solver<agent>, private generic_agent_solver<cell>
{
	// number of updates done so far; together with the seed it identifies the noise of a step
	index_t step_ = 0;

public:
	// forces holds agents_count * dims values, laid out like the positions
	void update(environment& e, std::span<const real_t> forces);
};

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
