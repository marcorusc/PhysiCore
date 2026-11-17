#include "micromechanics/cell_parameters.h"

#include <cmath>
#include <stdexcept>

namespace physicore::mechanics::micromechanics {

simulation_parameters scale_to_elements(const cell_parameters& cell, index_t elements_per_cell,
										const simulation_parameters& shape)
{
	if (elements_per_cell == 0)
	{
		throw std::invalid_argument("A cell needs at least one element");
	}
	if (cell.radius <= 0.0 || cell.stiffness <= 0.0 || cell.damping <= 0.0 || cell.packing_density <= 0.0)
	{
		throw std::invalid_argument("Cell radius, stiffness, damping and packing density must be positive");
	}
	if (shape.steepness <= 0.0 || shape.intra_well_depth <= 0.0)
	{
		throw std::invalid_argument("Steepness and intra-cellular well depth of the shape must be positive");
	}

	const auto count = static_cast<real_t>(elements_per_cell);
	const real_t inverse_cube_root = 1.0 / std::cbrt(count);

	const real_t pair_stiffness =
		cell.stiffness * inverse_cube_root * (1.0 - cell.scaling_correction * inverse_cube_root);
	if (pair_stiffness <= 0.0)
	{
		throw std::invalid_argument("Too few elements per cell for this scaling correction");
	}

	simulation_parameters params = shape;
	params.equilibrium_distance = 2.0 * cell.radius * std::cbrt(cell.packing_density / count);
	params.drag_coefficient = cell.damping / count;
	params.intra_well_depth = pair_stiffness * params.equilibrium_distance * params.equilibrium_distance
							  / (8.0 * shape.steepness * shape.steepness);
	params.inter_well_depth = params.intra_well_depth * (shape.inter_well_depth / shape.intra_well_depth);

	return params;
}

} // namespace physicore::mechanics::micromechanics
