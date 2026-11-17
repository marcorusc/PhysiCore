#include "micromechanics/cell_parameters.h"

#include <cmath>
#include <stdexcept>

namespace physicore::mechanics::micromechanics {

cell_type scale_to_elements(const cell_parameters& cell, index_t elements_per_cell, const simulation_parameters& params)
{
	if (elements_per_cell == 0)
	{
		throw std::invalid_argument("A cell needs at least one element");
	}
	if (cell.radius <= 0.0 || cell.stiffness <= 0.0 || cell.damping <= 0.0 || cell.packing_density <= 0.0)
	{
		throw std::invalid_argument("Cell radius, stiffness, damping and packing density must be positive");
	}
	if (params.steepness <= 0.0)
	{
		throw std::invalid_argument("The steepness of the potential must be positive");
	}
	if (cell.adhesion_ratio < 0.0)
	{
		throw std::invalid_argument("The adhesion ratio must not be negative");
	}

	const auto count = static_cast<real_t>(elements_per_cell);
	const real_t inverse_cube_root = 1.0 / std::cbrt(count);

	const real_t pair_stiffness =
		cell.stiffness * inverse_cube_root * (1.0 - cell.scaling_correction * inverse_cube_root);
	if (pair_stiffness <= 0.0)
	{
		throw std::invalid_argument("Too few elements per cell for this scaling correction");
	}

	cell_type type;
	type.elements_per_cell = elements_per_cell;
	type.equilibrium_distance = 2.0 * cell.radius * std::cbrt(cell.packing_density / count);
	type.drag_coefficient = cell.damping / count;
	type.intra_well_depth = pair_stiffness * type.equilibrium_distance * type.equilibrium_distance
							/ (8.0 * params.steepness * params.steepness);
	type.inter_well_depth = type.intra_well_depth * cell.adhesion_ratio;

	return type;
}

} // namespace physicore::mechanics::micromechanics
