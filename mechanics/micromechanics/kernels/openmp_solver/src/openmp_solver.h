#pragma once

#include <micromechanics/solver.h>

#include "cell_solver.h"
#include "force_solver.h"
#include "position_solver.h"

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

class openmp_solver : public solver
{
	force_solver forces;
	position_solver positions;
	cell_solver cells;

public:
	// Brings the cell-level quantities in line with the sub-agents
	void initialize(environment& e) override;

	// Advances the sub-agents and then updates the cell-level quantities once
	void solve(environment& e, index_t iterations) override;
};

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
