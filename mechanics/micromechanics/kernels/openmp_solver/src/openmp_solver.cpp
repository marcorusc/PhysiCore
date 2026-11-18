#include "openmp_solver.h"

#include <stdexcept>
#include <string>

#include <micromechanics/cell_type.h>

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

void openmp_solver::initialize(environment& e)
{
	check_cell_types(e);

	cells.update(e);
}

void openmp_solver::check_cell_types(environment& e)
{
	// building the table checks the values of every type
	const interaction_table table(e.cell_types, e.params);

	for (index_t c = 0; c < e.cells->size(); ++c)
	{
		if (e.cells->get_agent_at(c)->type_index() >= table.types_count())
		{
			throw std::invalid_argument("cell " + std::to_string(c) + " refers to a cell type that does not exist");
		}
	}
}

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
