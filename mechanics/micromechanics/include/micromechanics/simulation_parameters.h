#pragma once

#include <common/types.h>

#include "pair_potential.h"

namespace physicore::mechanics::micromechanics {

/**
 * @brief Parameters of the minimal sub-cellular element model.
 *
 * Elements of the same cell and elements of different cells interact through the same potential
 * shape and differ only in the well depth. The intra-cellular well must be the deeper one for a
 * cell to keep its integrity (Newman, Math. Biosci. Eng. 2005).
 *
 * steepness and cutoff_factor default to the values of SEM2 (Chattaraj et al. 2023, Table I).
 * The other defaults are unit values in arbitrary units, not a calibrated cell.
 */
struct simulation_parameters
{
	// Drag coefficient of one element: velocity = force / drag_coefficient
	real_t drag_coefficient = 1.0;

	// Distance at which two elements neither attract nor repel
	real_t equilibrium_distance = 1.0;

	// Dimensionless steepness (rho) of the potential; larger means stiffer and shorter ranged
	real_t steepness = 2.0;

	// Well depth for two elements of the same cell
	real_t intra_well_depth = 1.0;

	// Well depth for two elements of different cells
	real_t inter_well_depth = 0.5;

	// Interactions are ignored beyond cutoff_factor * equilibrium_distance
	real_t cutoff_factor = 2.5;

	pair_potential intra_potential() const
	{
		return { .well_depth = intra_well_depth,
				 .equilibrium_distance = equilibrium_distance,
				 .steepness = steepness,
				 .cutoff_distance = cutoff_factor * equilibrium_distance };
	}

	pair_potential inter_potential() const
	{
		return { .well_depth = inter_well_depth,
				 .equilibrium_distance = equilibrium_distance,
				 .steepness = steepness,
				 .cutoff_distance = cutoff_factor * equilibrium_distance };
	}
};

} // namespace physicore::mechanics::micromechanics
