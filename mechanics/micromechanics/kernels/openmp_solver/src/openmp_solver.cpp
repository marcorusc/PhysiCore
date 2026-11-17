#include "openmp_solver.h"

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

void openmp_solver::initialize(environment& /* e */)
{
	// The reference kernel keeps no state between timesteps.
}

void openmp_solver::solve(environment& e, index_t iterations)
{
	for (index_t i = 0; i < iterations; ++i)
	{
		forces.compute(e);

		positions.update(e, forces.forces());
	}
}

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
