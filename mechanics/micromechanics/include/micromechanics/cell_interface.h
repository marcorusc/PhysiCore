#pragma once

#include <cstdint>
#include <span>

#include <common/base_agent_interface.h>
#include <common/types.h>
#include <micromechanics/micromechanics_export.h>

namespace physicore::mechanics::micromechanics {

class MICROMECHANICS_EXPORT cell_interface : public virtual base_agent_interface
{
public:
	// Index of the cell's type in environment::cell_types
	virtual index_t& type_index() = 0;

	virtual std::uint8_t& is_movable() = 0;
	virtual std::uint8_t& is_motile() = 0;

	virtual real_t& migration_speed() = 0;
	virtual real_t& migration_bias() = 0;
	virtual std::span<real_t> migration_bias_direction() = 0;
	virtual real_t& persistence_time() = 0;
	virtual std::span<real_t> motility_direction() = 0;

	virtual std::span<real_t> velocity() = 0;
	virtual index_t agent_count() const = 0;
	// Pressure exerted on the cell by its surroundings: positive when compressed, negative when pulled
	// outwards. Set by the solver from the forces of sub-agents outside the cell.
	virtual real_t pressure() const = 0;

	virtual std::span<index_t> neighbors() = 0;
};

} // namespace physicore::mechanics::micromechanics
