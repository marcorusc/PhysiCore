#pragma once

#include <vector>

#include <common/generic_agent_solver.h>
#include <common/types.h>
#include <micromechanics/environment.h>

#include "neighbor_grid.h"

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

/**
 * @brief Pressure that the surroundings exert on every cell.
 *
 * For a cell c with centre X_c in d dimensions
 *
 *     pressure = -1 / (d * V_c) * sum over sub-agents i of c, sum over sub-agents j outside c,
 *                                 of (x_i - X_c) . F_ij
 *
 * with F_ij the force on i due to j. This is the virial of the outside forces. If the surroundings
 * press on the surface of the cell with a uniform pressure p, the sum is -p * d * V_c (divergence
 * theorem), so the expression returns p. A positive value means the cell is compressed, a negative one
 * that it is pulled outwards, as by adhesion to a neighbour that is being moved away.
 *
 * V_c is the nominal volume of the cell: its number of sub-agents times the volume one sub-agent takes
 * in a close packing at the equilibrium distance r_eq, which is r_eq^3 / sqrt(2) in three dimensions,
 * sqrt(3) / 2 * r_eq^2 in two and r_eq in one. For a cell made with scale_to_elements() and the default
 * packing density this is the volume of the sphere of the requested radius.
 *
 * Only forces from outside the cell count: from sub-agents of other cells and from sub-agents without a
 * cell. Forces inside the cell, the motility and the thermal noise do not. A cell without sub-agents, a
 * cell of one sub-agent (which sits on its own centre) and a cell without neighbours have zero
 * pressure. The value belongs to the positions at the time of the call; it is not an average over time.
 *
 * update() expects the cell positions to be current (see cell_solver). update_brute_force() gives the
 * same result up to rounding at quadratic cost and is the reference in the tests.
 */
class pressure_solver : private generic_agent_solver<agent>, private generic_agent_solver<cell>
{
	// contribution of every sub-agent to the virial of its cell
	std::vector<real_t> virials_;
	neighbor_grid grid_;

	void accumulate(environment& e);

public:
	void update(environment& e);

	void update_brute_force(environment& e);
};

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
