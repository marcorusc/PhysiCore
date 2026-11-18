#include "motility_solver.h"

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>

#include <micromechanics/counter_random.h>

#include "parallel_threshold.h"

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

namespace {
// The thermal noise uses the components 0 to 2 of the random counters, with the sub-agent index. The
// motility uses components from 8 on, with the cell index, so the two never share a number.
constexpr std::uint64_t turn_component = 8;
constexpr std::uint64_t direction_component = 9;

// below this length a vector is treated as having no direction
constexpr real_t no_length = 16.0 * std::numeric_limits<real_t>::epsilon();

real_t length_of(const std::array<real_t, 3>& vector)
{
	return std::sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
}

// Random unit vector, uniform on the circle (2D) or sphere (3D) and +-1 in 1D: independent normal
// numbers have no preferred direction, so normalising them gives a uniform one. Unused axes stay zero.
std::array<real_t, 3> random_direction(std::uint64_t seed, index_t step, index_t cell_index, index_t dims)
{
	std::array<real_t, 3> direction = { 0.0, 0.0, 0.0 };
	for (index_t d = 0; d < dims; ++d)
	{
		direction.at(d) = standard_normal(seed, step, cell_index, direction_component + d);
	}

	const real_t length = length_of(direction);
	if (length < no_length)
	{
		return { 1.0, 0.0, 0.0 };
	}
	for (real_t& component : direction)
	{
		component /= length;
	}
	return direction;
}
} // namespace

void motility_solver::validate(environment& e)
{
	const cell_data& cells = retrieve_agent_data(*e.cells);

	for (index_t c = 0; c < cells.cells_count; ++c)
	{
		if (cells.is_motile[c] == 0)
		{
			continue;
		}
		if (!(cells.persistence_times[c] > 0.0))
		{
			throw std::invalid_argument("motile cell " + std::to_string(c) + " needs a positive persistence time");
		}
		const real_t bias = cells.migration_biases[c];
		if (std::isnan(bias) || bias < 0.0 || bias > 1.0)
		{
			throw std::invalid_argument("motile cell " + std::to_string(c) + " needs a migration bias between 0 and 1");
		}
	}
}

void motility_solver::update(environment& e)
{
	validate(e);

	cell_data& cells = retrieve_agent_data(*e.cells);

	const index_t cells_count = cells.cells_count;
	const index_t dims = cells.dims;
	const real_t timestep = e.mechanics_timestep;
	const std::uint64_t seed = e.params.random_seed;
	const index_t step = step_;

	const std::uint8_t* is_motile = cells.is_motile.data();
	const real_t* persistence_times = cells.persistence_times.data();
	const real_t* biases = cells.migration_biases.data();
	const real_t* bias_directions = cells.migration_bias_directions.data();
	real_t* directions = cells.motility_directions.data();

#pragma omp parallel for if (cells_count >= parallel_threshold)
	for (index_t c = 0; c < cells_count; ++c)
	{
		if (is_motile[c] == 0)
		{
			continue;
		}

		std::array<real_t, 3> current = { 0.0, 0.0, 0.0 };
		for (index_t d = 0; d < dims; ++d)
		{
			current.at(d) = directions[c * dims + d];
		}
		const bool has_direction = length_of(current) >= no_length;

		// uniform in (0, 1], so the comparison is true with probability timestep / persistence_time
		const real_t uniform = to_unit_interval(random_bits(seed, step, c, turn_component, 0));
		const bool turns = uniform <= timestep / persistence_times[c];

		if (has_direction && !turns)
		{
			continue;
		}

		const std::array<real_t, 3> random = random_direction(seed, step, c, dims);

		std::array<real_t, 3> mixed = { 0.0, 0.0, 0.0 };
		for (index_t d = 0; d < dims; ++d)
		{
			mixed.at(d) = (1.0 - biases[c]) * random.at(d) + biases[c] * bias_directions[c * dims + d];
		}

		const real_t length = length_of(mixed);
		if (length < no_length)
		{
			continue;
		}
		for (index_t d = 0; d < dims; ++d)
		{
			directions[c * dims + d] = mixed.at(d) / length;
		}
	}

	++step_;
}

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
