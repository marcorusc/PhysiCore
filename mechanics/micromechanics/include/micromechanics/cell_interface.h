#pragma once

#include <cstdint>
#include <span>

#include <common/base_agent_interface.h>
#include <common/types.h>

namespace physicore::mechanics::micromechanics {

class cell_interface : public virtual base_agent_interface
{
public:
	// Index of the cell's type in environment::cell_types
	virtual index_t& type_index() = 0;

	virtual std::uint8_t& is_movable() = 0;

	virtual std::span<real_t> velocity() = 0;
	virtual index_t agent_count() const = 0;

	virtual std::span<index_t> neighbors() = 0;
};

} // namespace physicore::mechanics::micromechanics
