#pragma once

#include <memory>
#include <span>
#include <tuple>

#include <common/base_agent_data.h>
#include <common/base_agent_generic_storage.h>
#include <common/types.h>

#include "agent_data.h"
#include "agent_interface.h"

namespace physicore::mechanics::micromechanics {

#ifdef _MSC_VER
	#pragma warning(push)
	#pragma warning(disable : 4250) // inherits via dominance (intentional design)
#endif

class agent : public physicore::base_agent_generic_storage<base_agent_data>, public virtual agent_interface
{
protected:
	agent_data& data;

public:
	using DataType = agent_data;
	using InterfaceType = agent_interface;

	agent(index_t index, agent_data& data);
	agent(index_t index, std::tuple<std::unique_ptr<base_agent_data>, std::unique_ptr<agent_data>>& datas);

	index_t& cell_index() override;

	std::span<real_t> velocity() override;
};

#ifdef _MSC_VER
	#pragma warning(pop)
#endif

} // namespace physicore::mechanics::micromechanics
