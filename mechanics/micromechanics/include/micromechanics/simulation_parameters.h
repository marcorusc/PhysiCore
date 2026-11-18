#pragma once

#include <cstdint>

#include <common/types.h>

namespace physicore::mechanics::micromechanics {

/**
 * @brief Parameters of the sub-cellular element model that hold for the whole simulation.
 *
 * Everything that describes a kind of cell (size of its elements, how stiff they are, how strongly
 * cells stick to each other) is in cell_type. What is left here is the shape of the potential, which
 * all cell types share, the thermal noise and the walls of the domain.
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

	// Thermal energy kT of the noise. An element with the drag coefficient eta diffuses with
	// D = thermal_energy / eta (Einstein relation): a free element has the mean squared displacement
	// 2 * dims * D * time. One energy for all cell types means one temperature for all of them. Zero
	// switches the noise off and makes the dynamics deterministic. Compare it with the well depths to
	// judge whether bonds hold.
	real_t thermal_energy = 0.0;

	// Seed of the noise. The same seed gives the same simulation, bit for bit.
	std::uint64_t random_seed = 0;

	// Walls at the faces of the domain (the bounding box of the mesh) that push elements back in, see
	// domain_wall.h. Off by default: elements then move freely and the size of the domain plays no role.
	bool domain_walls = false;
};

} // namespace physicore::mechanics::micromechanics
