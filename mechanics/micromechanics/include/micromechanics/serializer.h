#pragma once

#include <memory>

#include <common/types.h>
#include <micromechanics/micromechanics_export.h>

namespace physicore::mechanics::micromechanics {

class environment;

class MICROMECHANICS_EXPORT serializer
{
public:
	serializer() = default;
	serializer(const serializer&) = delete;
	serializer& operator=(const serializer&) = delete;
	serializer(serializer&&) = delete;
	serializer& operator=(serializer&&) = delete;

	virtual void serialize(const environment& e, real_t current_time) = 0;

	virtual ~serializer() = default;
};

using serializer_ptr = std::unique_ptr<serializer>;

} // namespace physicore::mechanics::micromechanics
