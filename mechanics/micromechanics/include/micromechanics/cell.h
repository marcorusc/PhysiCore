#pragma once

#include <memory>
#include <tuple>

#include "cell_data.h"
#include "cell_interface.h"

namespace physicore::mechanics::micromechanics {

class cell : public virtual cell_interface
{
protected:
	cell_data& data;

public:
	using DataType = cell_data;
	using InterfaceType = cell_interface;

	cell(index_t index, cell_data& data);
	cell(index_t index, std::tuple<std::unique_ptr<cell_data>>& datas);

	std::span<real_t> position() override;

	index_t& type_index() override;
	std::uint8_t& is_movable() override;

	std::span<real_t> velocity() override;
	index_t agent_count() const override;

	std::span<index_t> neighbors() override;
};

} // namespace physicore::mechanics::micromechanics
