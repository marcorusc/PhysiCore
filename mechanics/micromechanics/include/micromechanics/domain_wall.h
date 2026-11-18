#pragma once

#include <algorithm>

#include <common/types.h>

#include "pair_potential.h"

namespace physicore::mechanics::micromechanics {

/**
 * Walls at the faces of the domain.
 *
 * An element at the distance h from a wall is pushed away from it as if it faced its own mirror image
 * behind the wall, at the distance 2h, through the pair potential of the wall. Only the repulsive part
 * acts: the force starts at h = equilibrium_distance / 2, where an element drawn as a sphere of that
 * radius touches the wall, and there is no attraction beyond it. The force is normal to the wall.
 *
 * The functions below treat one axis, with walls at `low` and `high`. An element outside the domain is
 * treated as sitting on the wall, where the force is zero: the repulsion of this potential falls back
 * to zero when two elements coincide. The solver therefore also puts back any element that leaves the
 * domain.
 */

/// Force along one axis on an element at `coordinate` between walls at `low` and `high`.
inline real_t wall_force(const pair_potential& wall, real_t coordinate, real_t low, real_t high)
{
	const real_t reach = 0.5 * wall.equilibrium_distance;

	real_t force = 0.0;

	const real_t from_low = std::max(coordinate - low, 0.0);
	if (from_low < reach)
	{
		// the image is at -from_low, so x_i - x_j = 2 * from_low; the coefficient is positive below r_eq
		force += force_coefficient(wall, 4.0 * from_low * from_low) * 2.0 * from_low;
	}

	const real_t from_high = std::max(high - coordinate, 0.0);
	if (from_high < reach)
	{
		force -= force_coefficient(wall, 4.0 * from_high * from_high) * 2.0 * from_high;
	}

	return force;
}

/// Potential energy of that force: wall_force is minus its derivative with respect to `coordinate`.
/// Zero where the wall does not act.
inline real_t wall_energy(const pair_potential& wall, real_t coordinate, real_t low, real_t high)
{
	const real_t reach = 0.5 * wall.equilibrium_distance;
	const real_t at_rest = potential_energy(wall, wall.equilibrium_distance * wall.equilibrium_distance);

	real_t energy = 0.0;
	for (const real_t distance : { std::max(coordinate - low, 0.0), std::max(high - coordinate, 0.0) })
	{
		if (distance < reach)
		{
			// half the pair energy at the distance 2h to the image, counted from its minimum
			energy += 0.5 * (potential_energy(wall, 4.0 * distance * distance) - at_rest);
		}
	}
	return energy;
}

} // namespace physicore::mechanics::micromechanics
