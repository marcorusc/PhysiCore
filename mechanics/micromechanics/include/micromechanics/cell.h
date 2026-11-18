#pragma once

#include <memory>
#include <tuple>

#include <micromechanics/micromechanics_export.h>

#include "cell_data.h"
#include "cell_interface.h"

namespace physicore::mechanics::micromechanics {

class MICROMECHANICS_EXPORT cell : public virtual cell_interface
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
	std::uint8_t& is_motile() override;

	real_t& migration_speed() override;
	real_t& migration_bias() override;
	std::span<real_t> migration_bias_direction() override;
	real_t& persistence_time() override;
	std::span<real_t> motility_direction() override;

	std::span<real_t> velocity() override;
	index_t agent_count() const override;
	real_t pressure() const override;

	std::span<index_t> neighbors() override;
};

} // namespace physicore::mechanics::micromechanics
