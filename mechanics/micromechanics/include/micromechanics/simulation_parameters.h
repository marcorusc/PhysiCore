#pragma once

#include <common/types.h>

namespace physicore::mechanics::micromechanics {

/**
 * @brief Parameters of the sub-cellular element model that hold for the whole simulation.
 *
 * Everything that describes a kind of cell (size of its elements, how stiff they are, how strongly
 * cells stick to each other) is in cell_type. What is left here is the shape of the potential, which
 * all cell types share.
 *
 * steepness defaults to the value used by Sandersius and Newman (2008) and SEM2 (Chattaraj et al.
 * 2023, Table I); cutoff_factor to the SEM2 value.
 */
struct simulation_parameters
{
	// Dimensionless steepness (rho) of the potential; larger means stiffer and shorter ranged
	real_t steepness = 2.0;

	// Two elements do not interact beyond cutoff_factor times their equilibrium distance
	real_t cutoff_factor = 2.5;
};

} // namespace physicore::mechanics::micromechanics
