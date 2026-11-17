#include "force_solver.h"

#include <array>

#include <micromechanics/domain_wall.h>
#include <micromechanics/pair_potential.h>

#include "parallel_threshold.h"

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

namespace {
// Everything needed to evaluate the force between two sub-agents
struct pair_forces
{
	const real_t* positions;
	const index_t* cell_indices;
	index_t dims;
	pair_potential intra;
	pair_potential inter;

	// Adds the force on sub-agent i due to sub-agent j to net_force
	void add(index_t i, index_t j, std::array<real_t, 3>& net_force) const
	{
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
};

pair_forces make_pair_forces(const agent_data& data, const simulation_parameters& params)
{
	return { .positions = data.base_data.positions.data(),
			 .cell_indices = data.cell_indices.data(),
			 .dims = data.base_data.dims,
			 .intra = params.intra_potential(),
			 .inter = params.inter_potential() };
}
} // namespace

void force_solver::compute(environment& e)
{
	const agent_data& data = retrieve_agent_data(*e.agents);
	const pair_forces pairs = make_pair_forces(data, e.params);

	const index_t agents_count = data.agents_count;
	const index_t dims = pairs.dims;

	// both potentials share the same cutoff distance
	grid_.build(pairs.positions, agents_count, dims, pairs.intra.cutoff_distance);

	forces_.resize(agents_count * dims);
	real_t* forces = forces_.data();
	const neighbor_grid& grid = grid_;

#pragma omp parallel for if (agents_count >= parallel_threshold)
	for (index_t i = 0; i < agents_count; ++i)
	{
		std::array<real_t, 3> net_force {};

		grid.visit_candidates(pairs.positions + i * dims, [&pairs, &net_force, i](index_t j) {
			if (i != j)
			{
				pairs.add(i, j, net_force);
			}
		});

		for (index_t d = 0; d < dims; ++d)
		{
			forces[i * dims + d] = net_force[d];
		}
	}

	add_wall_forces(e);
}

void force_solver::add_wall_forces(environment& e)
{
	if (!e.params.domain_walls)
	{
		return;
	}

	const agent_data& data = retrieve_agent_data(*e.agents);
	const index_t agents_count = data.agents_count;
	const index_t dims = data.base_data.dims;
	const pair_potential wall = e.params.wall_potential();
	const cartesian_mesh& mesh = e.mesh;

	const real_t* positions = data.base_data.positions.data();
	real_t* forces = forces_.data();

#pragma omp parallel for if (agents_count >= parallel_threshold)
	for (index_t i = 0; i < agents_count; ++i)
	{
		for (index_t d = 0; d < dims; ++d)
		{
			forces[i * dims + d] +=
				wall_force(wall, positions[i * dims + d], static_cast<real_t>(mesh.bounding_box_mins.at(d)),
						   static_cast<real_t>(mesh.bounding_box_maxs.at(d)));
		}
	}
}

void force_solver::compute_brute_force(environment& e)
{
	const agent_data& data = retrieve_agent_data(*e.agents);
	const pair_forces pairs = make_pair_forces(data, e.params);

	const index_t agents_count = data.agents_count;
	const index_t dims = pairs.dims;

	forces_.resize(agents_count * dims);
	real_t* forces = forces_.data();

#pragma omp parallel for if (agents_count >= parallel_threshold)
	for (index_t i = 0; i < agents_count; ++i)
	{
		std::array<real_t, 3> net_force {};

		for (index_t j = 0; j < agents_count; ++j)
		{
			if (i != j)
			{
				pairs.add(i, j, net_force);
			}
		}

		for (index_t d = 0; d < dims; ++d)
		{
			forces[i * dims + d] = net_force[d];
		}
	}

	add_wall_forces(e);
}

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
