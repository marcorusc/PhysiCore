#pragma once

#include <span>
#include <vector>

#include <common/generic_agent_solver.h>
#include <common/types.h>
#include <micromechanics/environment.h>

#include "neighbor_grid.h"

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

/**
 * @brief Net pair force on every sub-agent.
 *
 * compute() finds interaction partners through a uniform grid, so its cost grows linearly with the
 * number of sub-agents at fixed density. compute_brute_force() visits every pair and is the reference
 * the grid version is tested against.
 *
 * In both, the force on one sub-agent is summed in a fixed order by a single thread, which makes the
 * result independent of the number of threads. The two versions sum in different orders, so they agree
 * to rounding, not bitwise.
 */
class force_solver : private generic_agent_solver<agent>
{
	std::vector<real_t> forces_;
	neighbor_grid grid_;

public:
	// Recomputes the forces from the current positions
	void compute(environment& e);

	// Same result as compute() up to rounding, at quadratic cost
	void compute_brute_force(environment& e);

	// agents_count * dims values, laid out like the positions
	std::span<const real_t> forces() const { return forces_; }
};

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
