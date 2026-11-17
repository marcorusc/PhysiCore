#include "openmp_solver.h"

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

void openmp_solver::initialize(environment& e) { cells.update(e); }

void openmp_solver::solve(environment& e, index_t iterations)
{
	for (index_t i = 0; i < iterations; ++i)
	{
		forces.compute(e);

		positions.update(e, forces.forces());
	}

	cells.update(e);
}

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
