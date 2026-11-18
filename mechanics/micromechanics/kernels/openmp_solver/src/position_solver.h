#pragma once

#include <span>

#include <common/generic_agent_solver.h>
#include <common/types.h>
#include <micromechanics/environment.h>

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

/**
 * @brief Position update of the overdamped equation of motion (forward Euler scheme).
 *
 *     velocity = force / drag_coefficient
 *     position += velocity * timestep
 *
 * with drag_coefficient that of the sub-agent's cell type (Sandersius and Newman 2008, Eq. 1,
 * without the noise).
 *
 * Sub-agents of a cell that is not movable keep their position and get zero velocity.
 */
class position_solver : private generic_agent_solver<agent>, private generic_agent_solver<cell>
{
public:
	// forces holds agents_count * dims values, laid out like the positions
	void update(environment& e, std::span<const real_t> forces);
};

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
