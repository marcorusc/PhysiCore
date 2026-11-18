#include "micromechanics/solver_registry.h"

#include "openmp_solver/register_solver.h"

namespace physicore::mechanics::micromechanics {

solver_registry& solver_registry::instance()
{
	static solver_registry r;
	return r;
}

namespace {
struct attachment_point
{
	attachment_point() { kernels::openmp_solver::attach_to_registry(); }
};

const attachment_point ap;
} // namespace

} // namespace physicore::mechanics::micromechanics
