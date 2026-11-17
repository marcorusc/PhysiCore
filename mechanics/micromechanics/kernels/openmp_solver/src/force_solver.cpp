#include "force_solver.h"

#include <array>

#include <micromechanics/pair_potential.h>

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

void force_solver::compute(environment& e)
{
	const agent_data& data = retrieve_agent_data(*e.agents);

	const index_t agents_count = data.agents_count;
	const index_t dims = data.base_data.dims;

	const real_t* positions = data.base_data.positions.data();
	const index_t* cell_indices = data.cell_indices.data();

	const pair_potential intra = e.params.intra_potential();
	const pair_potential inter = e.params.inter_potential();

	forces_.resize(agents_count * dims);
	real_t* forces = forces_.data();

#pragma omp parallel for
	for (index_t i = 0; i < agents_count; ++i)
	{
		std::array<real_t, 3> net_force {};

		for (index_t j = 0; j < agents_count; ++j)
		{
			if (i == j)
				continue;

			std::array<real_t, 3> difference {};
			real_t squared_distance = 0.0;
			for (index_t d = 0; d < dims; ++d)
			{
				difference[d] = positions[i * dims + d] - positions[j * dims + d];
				squared_distance += difference[d] * difference[d];
			}

			const bool same_cell = cell_indices[i] == cell_indices[j] && cell_indices[i] != no_cell;
			const real_t coefficient = force_coefficient(same_cell ? intra : inter, squared_distance);

			for (index_t d = 0; d < dims; ++d)
			{
				net_force[d] += coefficient * difference[d];
			}
		}

		for (index_t d = 0; d < dims; ++d)
		{
			forces[i * dims + d] = net_force[d];
		}
	}
}

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
