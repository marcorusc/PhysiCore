#pragma once

#include <vector>

#include <common/base_agent_data.h>
#include <common/types.h>
#include <micromechanics/micromechanics_export.h>

namespace physicore::mechanics::micromechanics {

/**
 * @brief Structure-of-arrays storage for mechanical sub-agents.
 *
 * Positions are owned by the shared base_agent_data; this struct only refers to it.
 * The owning container calls add()/remove_at() on the base data first, so both stay in sync.
 */
struct MICROMECHANICS_EXPORT agent_data
{
	physicore::base_agent_data& base_data;

	index_t agents_count = 0;

	std::vector<index_t> cell_indices;

	std::vector<real_t> velocities; // agents_count * dims

	explicit agent_data(physicore::base_agent_data& base_data);

	void add();
	void remove_at(index_t position);
};

} // namespace physicore::mechanics::micromechanics
