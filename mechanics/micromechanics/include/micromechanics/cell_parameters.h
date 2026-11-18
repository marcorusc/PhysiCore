#pragma once

#include <numbers>

#include <common/types.h>

#include "cell_type.h"
#include "simulation_parameters.h"

namespace physicore::mechanics::micromechanics {

/**
 * @brief Measurable properties of one cell, independent of how many elements represent it.
 *
 * Following Sandersius and Newman, Phys. Biol. 5 (2008) 015002, section 2. Use one consistent set of
 * units for all values, for example micrometres, seconds and piconewtons: the paper's
 * stiffness 5e-3 N/m and damping 5e-3 N s/m are then 5000 pN/um and 5000 pN s/um.
 */
struct cell_parameters
{
	// Radius of the roughly spherical cell (R_cell)
	real_t radius = 0.0;

	// Stiffness scale (kappa_0); the elastic modulus of the cell is of the order stiffness / (2 * radius)
	real_t stiffness = 0.0;

	// Damping scale (eta_0); damping / stiffness is the time over which the cell relaxes
	real_t damping = 0.0;

	// Fraction of the cell volume filled by the elements (p_3); the default is close packing of spheres
	real_t packing_density = std::numbers::pi / (3.0 * std::numbers::sqrt2);

	// Leading correction to the stiffness scaling (lambda); 0.75 is the value used in the paper
	real_t scaling_correction = 0.75;

	// Depth of the well between two cells of this kind, as a fraction of the well inside the cell.
	// Below 1 the cells adhere and stay distinct; at 1 nothing tells two cells apart and they merge.
	real_t adhesion_ratio = 0.5;
};

/**
 * @brief The cell type whose elements make a cell of elements_per_cell elements behave like cell.
 *
 * With N = elements_per_cell (the paper recommends N of about 100 or more):
 *
 *     equilibrium_distance = 2 * radius * (packing_density / N)^(1/3)                     (Eq. 8)
 *     pair stiffness kappa = stiffness * N^(-1/3) * (1 - scaling_correction * N^(-1/3))   (Eq. 11)
 *     drag_coefficient     = damping / N                                                  (Eq. 12)
 *     intra_well_depth     = kappa * equilibrium_distance^2 / (8 * steepness^2)           (Eq. 5)
 *     inter_well_depth     = adhesion_ratio * intra_well_depth
 *
 * The steepness is taken from params, which must be the parameters the simulation runs with. The
 * relations are those of the paper and hold in three dimensions only. Name and affinities
 * of the returned type are left at their defaults.
 *
 * @throws std::invalid_argument if a value is not positive or the well depth would not be positive
 */
cell_type scale_to_elements(const cell_parameters& cell, index_t elements_per_cell,
							const simulation_parameters& params = {});

} // namespace physicore::mechanics::micromechanics
