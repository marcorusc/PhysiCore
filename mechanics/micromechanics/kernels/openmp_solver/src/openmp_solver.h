#pragma once

#include <micromechanics/solver.h>

#include "cell_solver.h"
#include "force_solver.h"
#include "motility_solver.h"
#include "position_solver.h"

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

class openmp_solver : public solver
{
	motility_solver motility;
	force_solver forces;
	position_solver positions;
	cell_solver cells;

	// Throws std::invalid_argument if a cell type has an invalid value or a cell refers to a type
	// that does not exist
	static void check_cell_types(environment& e);

public:
	// Checks the cell types and the motility parameters, and brings the cell-level quantities in line
	// with the sub-agents
	void initialize(environment& e) override;

	// Every iteration chooses the motility directions, computes the forces and moves the sub-agents.
	// The cell-level quantities are updated once at the end.
	void solve(environment& e, index_t iterations) override;
};

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
