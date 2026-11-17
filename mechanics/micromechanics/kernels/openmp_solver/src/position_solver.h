#pragma once

#include <span>

#include <common/generic_agent_solver.h>
#include <common/types.h>
#include <micromechanics/environment.h>

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

/**
 * @brief Overdamped position update with forward Euler.
 *
 *     velocity = force / drag_coefficient
 *     position += velocity * mechanics_timestep
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
