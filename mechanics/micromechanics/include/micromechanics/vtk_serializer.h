#pragma once

#include <string_view>

#include <micromechanics/micromechanics_export.h>

#include "serializer.h"

namespace physicore::mechanics::micromechanics {

/**
 * @brief Creates a serializer that writes the sub-agents as VTK files readable by ParaView.
 *
 * Every call to serialize() writes <output_dir>/vtk_micromechanics_agents/micromechanics_agents_NNNNNN.vtu
 * and updates the collection file <output_dir>/micromechanics_agents.pvd, which lists all written files
 * with their times. Each sub-agent is one point carrying:
 *
 * - cell_index: index of its cell, -1 if it has none (color by this to tell cells apart)
 * - cell_type: index of the type of its cell, -1 if it has no cell
 * - velocity: three components, padded with zeros below three dimensions
 */
MICROMECHANICS_EXPORT serializer_ptr make_vtk_agents_serializer(std::string_view output_dir);

} // namespace physicore::mechanics::micromechanics
