#include "micromechanics/solver_registry.h"

namespace physicore::mechanics::micromechanics {

solver_registry& solver_registry::instance()
{
	static solver_registry r;
	return r;
}

} // namespace physicore::mechanics::micromechanics
