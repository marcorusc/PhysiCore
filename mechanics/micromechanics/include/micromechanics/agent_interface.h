#pragma once

#include <limits>
#include <span>

#include <common/base_agent_interface.h>
#include <common/types.h>
#include <micromechanics/micromechanics_export.h>

namespace physicore::mechanics::micromechanics {

/// Value of agent_interface::cell_index() for a sub-agent that is not assigned to any cell.
inline constexpr index_t no_cell = std::numeric_limits<index_t>::max();

/**
 * @brief A mechanical sub-agent: one element of a sub-cellular element model cell.
 *
 * Several sub-agents together represent one biological cell. Properties shared by all
 * sub-agents of a cell (motility, ...) live on the cell, see cell_interface.
 */
class MICROMECHANICS_EXPORT agent_interface : public virtual base_agent_interface
{
public:
	/// Index of the owning cell in the cell container, or no_cell.
	virtual index_t& cell_index() = 0;

	virtual std::span<real_t> velocity() = 0;
};

} // namespace physicore::mechanics::micromechanics
