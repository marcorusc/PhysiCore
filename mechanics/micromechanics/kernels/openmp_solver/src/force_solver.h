#pragma once

#include <span>
#include <vector>

#include <common/generic_agent_solver.h>
#include <common/types.h>
#include <micromechanics/environment.h>

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

/**
 * @brief Net pair force on every sub-agent.
 *
 * Reference implementation: every sub-agent visits every other one, so the cost is quadratic in the
 * number of sub-agents. Each force is summed in a fixed order by a single thread, which makes the
 * result independent of the number of threads.
 */
class force_solver : private generic_agent_solver<agent>
{
	std::vector<real_t> forces_;

public:
	// Recomputes the forces from the current positions
	void compute(environment& e);

	// agents_count * dims values, laid out like the positions
	std::span<const real_t> forces() const { return forces_; }
};

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
