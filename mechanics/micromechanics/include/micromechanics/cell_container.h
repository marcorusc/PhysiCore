#pragma once

#include <memory>

#include <common/generic_agent_container.h>

#include "cell.h"

namespace physicore::mechanics::micromechanics {

using cell_container = physicore::generic_agent_and_data_container<cell>;
using cell_container_interface = physicore::generic_agent_interface_container<cell_interface>;
using cell_container_ptr = std::shared_ptr<cell_container_interface>;

} // namespace physicore::mechanics::micromechanics
