#pragma once

#include <algorithm>
#include <array>
#include <initializer_list>
#include <random>
#include <vector>

#include <micromechanics/cell_type.h>
#include <micromechanics/environment.h>
#include <micromechanics/pair_potential.h>

#include "force_solver.h"

namespace physicore::mechanics::micromechanics::kernels::openmp_solver::tests {

inline environment make_environment(index_t dims, real_t timestep = 0.001)
{
	return { { dims, { -100, -100, -100 }, { 100, 100, 100 }, { 20, 20, 20 } }, timestep };
}

// Net force on every sub-agent, agents_count * dims values
inline std::vector<real_t> forces_of(environment& env)
{
	force_solver forces;
	forces.compute(env);
	return { forces.forces().begin(), forces.forces().end() };
}

// A cell of one sub-agent of the given type on the x axis
inline cell_interface* add_point_cell(environment& env, index_t type_index, real_t x)
{
	auto* cell = env.create_cell(type_index);
	env.create_agent(env.cells->size() - 1)->position()[0] = x;
	return cell;
}

inline real_t x_of(environment& env, index_t agent) { return env.agents->get_agent_at(agent)->position()[0]; }

// A cell type whose elements are twice as large as those of the default type
inline cell_type large_type()
{
	cell_type type;
	type.name = "large";
	type.equilibrium_distance = 2.0;
	return type;
}

inline agent_interface* add_agent(environment& env, index_t cell_index, std::initializer_list<real_t> position)
{
	auto* created = env.create_agent(cell_index);
	std::ranges::copy(position, created->position().begin());
	return created;
}

// Copy of all sub-agent positions, agents_count * dims values
inline std::vector<real_t> positions_of(environment& env)
{
	std::vector<real_t> positions;
	for (index_t i = 0; i < env.agents->size(); ++i)
	{
		const auto position = env.agents->get_agent_at(i)->position();
		positions.insert(positions.end(), position.begin(), position.end());
	}
	return positions;
}

inline std::array<real_t, 3> centre_of_mass(environment& env)
{
	std::array<real_t, 3> centre {};
	for (index_t i = 0; i < env.agents->size(); ++i)
	{
		const auto position = env.agents->get_agent_at(i)->position();
		for (std::size_t d = 0; d < position.size(); ++d)
		{
			centre.at(d) += position[d] / static_cast<real_t>(env.agents->size());
		}
	}
	return centre;
}

inline real_t squared_distance(environment& env, index_t i, index_t j)
{
	const auto a = env.agents->get_agent_at(i)->position();
	const auto b = env.agents->get_agent_at(j)->position();

	real_t result = 0.0;
	for (std::size_t d = 0; d < a.size(); ++d)
	{
		result += (a[d] - b[d]) * (a[d] - b[d]);
	}
	return result;
}

// Total potential energy, each pair counted once
inline real_t total_energy(environment& env)
{
	const pair_potential intra = intra_potential(env.cell_types[0], env.params);
	const pair_potential inter = inter_potential(env.cell_types[0], env.params);

	real_t energy = 0.0;
	for (index_t i = 0; i < env.agents->size(); ++i)
	{
		for (index_t j = i + 1; j < env.agents->size(); ++j)
		{
			const index_t cell_i = env.agents->get_agent_at(i)->cell_index();
			const index_t cell_j = env.agents->get_agent_at(j)->cell_index();
			const bool same_cell = cell_i == cell_j && cell_i != no_cell;

			energy += potential_energy(same_cell ? intra : inter, squared_distance(env, i, j));
		}
	}
	return energy;
}

// One cell of 27 sub-agents on a 3 x 3 x 3 lattice with the given spacing, each displaced by up to
// +-jitter along every axis. The lattice starts at x = x_offset. The seed is fixed, so the configuration
// is the same on every run.
inline void add_jittered_lattice_cell(environment& env, real_t spacing, real_t jitter, real_t x_offset = 0.0,
									  unsigned seed = 42)
{
	std::mt19937 generator(seed);
	std::uniform_real_distribution<real_t> displacement(-jitter, jitter);

	env.create_cell();
	const index_t cell_index = env.cells->size() - 1;

	for (int x = 0; x < 3; ++x)
	{
		for (int y = 0; y < 3; ++y)
		{
			for (int z = 0; z < 3; ++z)
			{
				add_agent(env, cell_index,
						  { x_offset + x * spacing + displacement(generator), y * spacing + displacement(generator),
							z * spacing + displacement(generator) });
			}
		}
	}
}

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver::tests
