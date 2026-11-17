#include "openmp_solver.h"

#include <stdexcept>
#include <string>

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

void openmp_solver::initialize(environment& e)
{
	motility.validate(e);
	check_domain(e);

	cells.update(e);
	pressures.update(e);
}

void openmp_solver::check_domain(environment& e)
{
	if (!e.params.domain_walls)
	{
		return;
	}

	const index_t dims = e.mesh.dims;
	for (index_t i = 0; i < e.agents->size(); ++i)
	{
		const auto position = e.agents->get_agent_at(i)->position();
		for (index_t d = 0; d < dims; ++d)
		{
			if (position[d] < static_cast<real_t>(e.mesh.bounding_box_mins.at(d))
				|| position[d] > static_cast<real_t>(e.mesh.bounding_box_maxs.at(d)))
			{
				throw std::invalid_argument("sub-agent " + std::to_string(i)
											+ " is outside the domain, which has walls");
			}
		}
	}
}

void openmp_solver::solve(environment& e, index_t iterations)
{
	for (index_t i = 0; i < iterations; ++i)
	{
		motility.update(e);

		forces.compute(e);

		positions.update(e, forces.forces());
	}

	cells.update(e);
	pressures.update(e);
}

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
