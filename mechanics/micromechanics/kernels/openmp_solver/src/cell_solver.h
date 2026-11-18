#pragma once

#include <common/generic_agent_solver.h>
#include <micromechanics/environment.h>

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

/**
 * @brief Derives the cell-level quantities from the sub-agents.
 *
 * For every cell: the number of sub-agents, the position (mean sub-agent position, i.e. the centre of
 * mass of identical elements) and the velocity (mean sub-agent velocity). A cell without sub-agents
 * keeps its position and velocity.
 */
class cell_solver : private generic_agent_solver<agent>, private generic_agent_solver<cell>
{
public:
	void update(environment& e);
};

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
