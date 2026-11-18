#pragma once

#include <span>
#include <string>
#include <vector>

#include <common/types.h>

#include "pair_potential.h"
#include "simulation_parameters.h"

namespace physicore::mechanics::micromechanics {

/**
 * @brief A kind of cell: what its elements are like and how its cells stick to others.
 *
 * Every cell refers to one cell type. The defaults are unit values in arbitrary units, not a calibrated
 * cell; scale_to_elements() (cell_parameters.h) computes the element values from measurable properties
 * of a cell.
 *
 * Two elements interact in one of three ways, see interaction_table:
 *
 * - elements of the same cell: through the intra-cellular well of the cell's type;
 * - elements of two cells of the same type: through the inter-cellular well of that type;
 * - elements of two cells of different types: through a well mixed from the two types.
 *
 * Between cells, the attraction is weighted by adhesion affinities, as in PhysiCell.
 */
struct cell_type
{
	std::string name = "default";

	// Drag coefficient of one element: velocity = force / drag_coefficient
	real_t drag_coefficient = 1.0;

	// Distance at which two elements of this type neither attract nor repel; an element can be
	// pictured as a sphere of half that diameter
	real_t equilibrium_distance = 1.0;

	// Well depth for two elements of the same cell
	real_t intra_well_depth = 1.0;

	// Well depth for two elements of different cells of this type. It must be the shallower of the
	// two for cells to stay distinct (Newman, Math. Biosci. Eng. 2005).
	real_t inter_well_depth = 0.5;

	// adhesion_affinities[s], between 0 and 1, weights the attraction to cells of type s; the entry
	// for this type itself weights the attraction between two cells of this type. Types without an
	// entry have the affinity 1. An affinity changes the attraction only: two cells repel each other
	// in full whatever their affinity, so cells with the affinity 0 still cannot overlap.
	std::vector<real_t> adhesion_affinities;

	// Number of elements of a cell created without giving a number
	index_t elements_per_cell = 1;

	// Affinity to type other_type, 1 if there is no entry for it
	real_t affinity_to(index_t other_type) const
	{
		return other_type < adhesion_affinities.size() ? adhesion_affinities[other_type] : 1.0;
	}
};

/// Potential between two elements of one cell of this type
inline pair_potential intra_potential(const cell_type& type, const simulation_parameters& params)
{
	return { .well_depth = type.intra_well_depth,
			 .equilibrium_distance = type.equilibrium_distance,
			 .steepness = params.steepness,
			 .cutoff_distance = params.cutoff_factor * type.equilibrium_distance };
}

/// Potential between two elements of different cells of this type, before the affinity is applied
inline pair_potential inter_potential(const cell_type& type, const simulation_parameters& params)
{
	return { .well_depth = type.inter_well_depth,
			 .equilibrium_distance = type.equilibrium_distance,
			 .steepness = params.steepness,
			 .cutoff_distance = params.cutoff_factor * type.equilibrium_distance };
}

/**
 * @brief A pair potential whose attractive side is weighted.
 *
 * Closer than the equilibrium distance the two elements repel as the potential says. Further apart,
 * the attraction is multiplied by attraction_scale, between 0 and 1. The force is zero at the
 * equilibrium distance in both cases, so it stays continuous there and the equilibrium distance does
 * not depend on the scale. With the scale 1 this is the plain potential; with 0 it only repels.
 */
struct pair_interaction
{
	pair_potential potential;
	real_t attraction_scale = 1.0;
};

/// Force coefficient, with the sign convention of force_coefficient(const pair_potential&, real_t)
inline real_t force_coefficient(const pair_interaction& interaction, real_t squared_distance)
{
	const real_t coefficient = force_coefficient(interaction.potential, squared_distance);
	return coefficient < 0.0 ? interaction.attraction_scale * coefficient : coefficient;
}

/// Potential energy of which that force is minus the gradient; continuous at the equilibrium distance
inline real_t potential_energy(const pair_interaction& interaction, real_t squared_distance)
{
	const pair_potential& potential = interaction.potential;
	const real_t energy = potential_energy(potential, squared_distance);

	if (squared_distance >= potential.equilibrium_distance * potential.equilibrium_distance)
	{
		return interaction.attraction_scale * energy;
	}
	// the repulsive side keeps its shape and is shifted to meet the weighted well at its bottom
	return energy + (1.0 - interaction.attraction_scale) * potential.well_depth;
}

/**
 * @brief The interaction of every pair of cell types, worked out once from the types.
 *
 * For two elements of different cells, of the types t and s:
 *
 *     equilibrium distance = (r_t + r_s) / 2          the sum of the two element radii
 *     well depth           = sqrt(c_t * c_s)          c: inter_well_depth of the type
 *     attraction scale     = sqrt(a_ts * a_st)        a_ts: affinity of type t to type s
 *
 * The geometric means are those of PhysiCell. They make the interaction the same seen from either
 * element, whatever the two affinity entries are, so action equals reaction. For two cells of the same
 * type the rule gives the inter-cellular well of that type, weighted by its affinity to itself.
 */
class interaction_table
{
	index_t types_count_ = 0;
	std::vector<pair_interaction> intra_; // one per type
	std::vector<pair_interaction> inter_; // types_count * types_count
	real_t largest_cutoff_ = 0.0;

public:
	/// @throws std::invalid_argument if there is no type, a type has a value that is not positive
	///         (the inter-cellular well depth may be zero) or an affinity outside [0, 1]
	interaction_table(std::span<const cell_type> types, const simulation_parameters& params);

	index_t types_count() const { return types_count_; }

	// two elements of the same cell
	const pair_interaction& intra(index_t type) const { return intra_[type]; }

	// two elements of different cells
	const pair_interaction& inter(index_t type, index_t other_type) const
	{
		return inter_[type * types_count_ + other_type];
	}

	// no two elements interact beyond this distance
	real_t largest_cutoff() const { return largest_cutoff_; }
};

} // namespace physicore::mechanics::micromechanics
