#pragma once

#include <array>
#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <common/mesh.h>
#include <common/types.h>
#include <micromechanics/micromechanics_export.h>

#include "cell_type.h"
#include "environment.h"
#include "simulation_parameters.h"

namespace physicore::mechanics::micromechanics {

/**
 * @brief Largest timestep at which forward Euler is stable for elements of this type.
 *
 * An element held by springs of stiffness kappa relaxes with rates up to lambda / drag_coefficient,
 * where lambda is the largest eigenvalue of the stiffness matrix of the packing. Forward Euler is
 * stable for timestep * lambda / drag_coefficient < 2. For a close packing with springs between nearest
 * neighbours lambda is 4 kappa on a line, 6 kappa in the hexagonal packing in the plane and 8 kappa in
 * the face-centred cubic packing in space, so
 *
 *     timestep < drag_coefficient / (2 kappa), / (3 kappa), / (4 kappa)      in 1, 2, 3 dimensions
 *
 * kappa = 8 * well_depth * steepness^2 / equilibrium_distance^2 is the stiffness at the bottom of the
 * deeper of the two wells of the type (Sandersius and Newman 2008, Eq. 5). Beyond the limit the
 * elements of a cell oscillate with growing amplitude and the cell flies apart. Stay well below it:
 * the examples use a fifth of it or less.
 */
MICROMECHANICS_EXPORT real_t largest_stable_timestep(const cell_type& type, const simulation_parameters& params,
													 index_t dims);

/**
 * @brief Puts an environment together step by step and checks it before it is used.
 *
 * build() returns an environment with its cell types, parameters and solver in place and the solver
 * initialised. It refuses values that cannot work: no domain, a cell type with a value that makes no
 * sense, or a timestep beyond largest_stable_timestep() of any cell type.
 */
class MICROMECHANICS_EXPORT environment_builder
{
	real_t timestep = 0.0;
	std::optional<cartesian_mesh> mesh;
	simulation_parameters params;
	std::vector<cell_type> cell_types;
	std::string solver_name = "openmp_solver";
	std::string vtk_output;

public:
	void set_time_step(real_t time_step);

	// the domain: its bounding box, and the voxel size of the mesh shared with other modules
	void resize(index_t dims, std::array<sindex_t, 3> bounding_box_mins, std::array<sindex_t, 3> bounding_box_maxs,
				std::array<index_t, 3> voxel_shape);

	void set_parameters(const simulation_parameters& parameters);
	void set_domain_walls(bool value);

	// Adds a cell type and returns its index. Without any, the environment has the default type.
	index_t add_cell_type(cell_type type);

	void set_solver(std::string name);

	// Writes the sub-agents as VTK files into this folder whenever the state is serialized
	void set_vtk_output(std::string output_dir);

	/// @throws std::invalid_argument if the domain or the timestep is missing, a cell type is invalid,
	///         the timestep is not stable or the solver is not registered
	std::unique_ptr<environment> build();
};

} // namespace physicore::mechanics::micromechanics
