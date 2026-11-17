#pragma once

#include <cstdint>

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
 * steepness defaults to the value used by Sandersius and Newman (2008) and SEM2 (Chattaraj et al.
 * 2023, Table I); cutoff_factor to the SEM2 value.
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

	// Strength D of the thermal noise: a free element diffuses with mean squared displacement
	// 2 * dims * D * time. Zero switches the noise off and makes the dynamics deterministic.
	// D * drag_coefficient plays the role of the thermal energy kT.
	real_t diffusion_coefficient = 0.0;

	// Seed of the noise. The same seed gives the same simulation, bit for bit.
	std::uint64_t random_seed = 0;

	// Walls at the faces of the domain (the bounding box of the mesh) that push elements back in, see
	// domain_wall.h. Off by default: elements then move freely and the size of the domain plays no role.
	bool domain_walls = false;

	// Well depth of the potential between an element and its mirror image behind a wall. The default
	// is that of two elements of different cells.
	real_t wall_well_depth = 0.5;

	pair_potential intra_potential() const
	{
		return { .well_depth = intra_well_depth,
				 .equilibrium_distance = equilibrium_distance,
				 .steepness = steepness,
				 .cutoff_distance = cutoff_factor * equilibrium_distance };
	}

	pair_potential wall_potential() const
	{
		return { .well_depth = wall_well_depth,
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
