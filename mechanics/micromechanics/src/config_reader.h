#pragma once

#include <array>
#include <filesystem>
#include <optional>
#include <vector>

#include <common/types.h>

#include "micromechanics/cell_type.h"
#include "micromechanics/simulation_parameters.h"

namespace physicore::mechanics::micromechanics {

/// What a PhysiCell settings file says about a sub-cellular element simulation
struct simulation_config
{
	// <domain>
	index_t dims = 3;
	std::array<sindex_t, 3> bounding_box_mins {};
	std::array<sindex_t, 3> bounding_box_maxs {};
	std::array<index_t, 3> voxel_shape {};

	// <overall><dt_mechanics>
	real_t timestep = 0.0;

	// <micromechanics>, and <options><virtual_wall_at_domain_edge> for the walls
	simulation_parameters params;

	// one per <cell_definition>, at the index given by its ID
	std::vector<cell_type> cell_types;

	// <initial_conditions><cell_positions>, if enabled
	std::optional<std::filesystem::path> cell_positions_file;
};

/**
 * @brief Reads the settings of a sub-cellular element simulation from a PhysiCell settings file.
 *
 * The file is the PhysiCell_settings.xml that PhysiCell itself reads. Fields it already has are used
 * where they mean the same here:
 *
 *     <domain>                                   size of the domain, 2D or 3D
 *     <overall><dt_mechanics>                    timestep
 *     <options><virtual_wall_at_domain_edge>     walls at the faces of the domain (off if absent)
 *     <cell_definition name ID>                  a cell type
 *       <phenotype><volume><total>               its volume, hence the radius of the cell
 *       <phenotype><mechanics>
 *         <cell_adhesion_affinities>             its affinities, by the name of the other type (1 if absent)
 *       <phenotype><motility>                    speed, persistence_time, migration_bias, options/enabled
 *     <initial_conditions><cell_positions>       a file with one cell per line
 *
 * Two blocks are specific to this model. In every cell definition, under <phenotype><mechanics>:
 *
 *     <sub_cellular_elements>
 *       <elements>128</elements>                 number of elements of a cell of this type
 *       <stiffness>5000</stiffness>              of the whole cell, see cell_parameters
 *       <damping>5000</damping>
 *       <adhesion_ratio>0.5</adhesion_ratio>     optional, 0.5 if absent
 *     </sub_cellular_elements>
 *
 * from which scale_to_elements() works out the elements; that scaling holds in three dimensions only.
 * Alternatively the elements are given directly, which works in any dimension:
 *
 *     <sub_cellular_elements>
 *       <elements>19</elements>
 *       <element>
 *         <equilibrium_distance>1</equilibrium_distance>
 *         <intra_well_depth>1</intra_well_depth>
 *         <inter_well_depth>0.5</inter_well_depth>
 *         <drag_coefficient>1</drag_coefficient>
 *       </element>
 *     </sub_cellular_elements>
 *
 * And at the top level, optional, for what all cell types share (defaults of simulation_parameters):
 *
 *     <micromechanics>
 *       <steepness>2</steepness>
 *       <cutoff_factor>2.5</cutoff_factor>
 *       <thermal_energy>0</thermal_energy>
 *       <random_seed>0</random_seed>
 *     </micromechanics>
 *
 * Everything else in the file is ignored, among it the cycle, death and volume dynamics of the
 * phenotype, PhysiCell's own adhesion and repulsion strengths and chemotaxis. Units are not
 * interpreted: all values must be in one consistent set of units.
 *
 * @throws std::runtime_error if the file cannot be read, a required tag is missing, a value makes no
 *         sense, an affinity names an unknown cell type, IDs are not 0, 1, 2, ... or a cell type gives
 *         cell-level values in a two-dimensional simulation
 */
simulation_config parse_simulation_config(const std::filesystem::path& config_file);

} // namespace physicore::mechanics::micromechanics
