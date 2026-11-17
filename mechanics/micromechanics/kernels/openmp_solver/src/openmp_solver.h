#pragma once

#include <micromechanics/solver.h>

#include "cell_solver.h"
#include "force_solver.h"
#include "motility_solver.h"
#include "position_solver.h"
#include "pressure_solver.h"

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

class openmp_solver : public solver
{
	motility_solver motility;
	force_solver forces;
	position_solver positions;
	cell_solver cells;
	pressure_solver pressures;

	// Throws std::invalid_argument if a cell type has an invalid value or a cell refers to a type
	// that does not exist
	static void check_cell_types(environment& e);

	// Throws std::invalid_argument if the domain has walls and a sub-agent is outside it
	static void check_domain(environment& e);

public:
	// Number of times a sub-agent had to be put back on a wall of the domain (see position_solver)
	index_t wall_stops() const { return positions.wall_stops(); }

	// Checks the cell types and the motility parameters, and brings the cell-level quantities in line
	// with the sub-agents
	void initialize(environment& e) override;

	// Every iteration chooses the motility directions, computes the forces and moves the sub-agents.
	// The cell-level quantities, the pressure among them, are updated once at the end.
	void solve(environment& e, index_t iterations) override;
};

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
