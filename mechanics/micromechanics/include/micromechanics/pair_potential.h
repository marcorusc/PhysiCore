#pragma once

#include <cmath>

#include <common/types.h>

namespace physicore::mechanics::micromechanics {

/**
 * @brief Short-range element-element potential of the sub-cellular element model.
 *
 * With t(r) = exp(steepness * (1 - r^2 / equilibrium_distance^2)):
 *
 *     V(r) = well_depth * (t^2 - 2 t)        for r < cutoff_distance, 0 otherwise
 *
 * (Sandersius and Newman, Phys. Biol. 5 (2008) 015002, Eq. 4). The minimum is at
 * r = equilibrium_distance with V = -well_depth, and the stiffness there is
 * V'' = 8 * well_depth * steepness^2 / equilibrium_distance^2 (Eq. 5 of the same paper).
 *
 * The potential is not shifted at the cutoff; with the default parameters the jump is below
 * 1e-4 * well_depth.
 */
struct pair_potential
{
	real_t well_depth;
	real_t equilibrium_distance;
	real_t steepness;
	real_t cutoff_distance;
};

/// Potential energy V of two elements whose squared distance is squared_distance.
inline real_t potential_energy(const pair_potential& potential, real_t squared_distance)
{
	if (squared_distance >= potential.cutoff_distance * potential.cutoff_distance)
		return 0.0;

	const real_t scaled = squared_distance / (potential.equilibrium_distance * potential.equilibrium_distance);
	const real_t t = std::exp(potential.steepness * (1.0 - scaled));

	return potential.well_depth * t * (t - 2.0);
}

/**
 * @brief Force coefficient c of two elements whose squared distance is squared_distance.
 *
 * Sign convention, used by every backend: the force on element i due to element j is
 *
 *     F_i = c * (x_i - x_j)          and          F_j = -F_i
 *
 * so c > 0 pushes the two elements apart (repulsion, r < equilibrium_distance) and c < 0 pulls them
 * together (adhesion). c equals -V'(r) / r; it needs neither a square root nor a division by r and is
 * finite at r = 0.
 */
inline real_t force_coefficient(const pair_potential& potential, real_t squared_distance)
{
	if (squared_distance >= potential.cutoff_distance * potential.cutoff_distance)
		return 0.0;

	const real_t inverse_squared_equilibrium = 1.0 / (potential.equilibrium_distance * potential.equilibrium_distance);
	const real_t t = std::exp(potential.steepness * (1.0 - squared_distance * inverse_squared_equilibrium));

	return 4.0 * potential.well_depth * potential.steepness * inverse_squared_equilibrium * t * (t - 1.0);
}

} // namespace physicore::mechanics::micromechanics
