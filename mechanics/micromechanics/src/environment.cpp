#include "micromechanics/environment.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

#include <common/base_agent_data.h>

#include "micromechanics/cell_placement.h"
#include "micromechanics/counter_random.h"

namespace physicore::mechanics::micromechanics {

environment::environment(const cartesian_mesh& mesh, real_t timestep) : mechanics_timestep(timestep), mesh(mesh)
{
	cells = std::make_shared<cell_container>(std::make_unique<cell_data>(mesh.dims));

	auto base_data = std::make_unique<physicore::base_agent_data>(mesh.dims);
	auto data = std::make_unique<agent_data>(*base_data);
	agents = std::make_shared<agent_container>(std::move(base_data), std::move(data));
}

void environment::run_single_timestep()
{
	if (!solver)
	{
		throw std::runtime_error("Micromechanics solver is not set");
	}

	solver->solve(*this, 1);
}

void environment::serialize_state(real_t current_time)
{
	if (serializer)
	{
		serializer->serialize(*this, current_time);
	}
}

cell_interface* environment::create_cell() { return cells->create(); }

cell_interface* environment::create_cell(std::span<const real_t> centre, index_t elements_count)
{
	const index_t dims = mesh.dims;
	if (centre.size() != dims)
	{
		throw std::invalid_argument("The centre of a cell needs one value per dimension");
	}

	const std::vector<real_t> offsets = packed_ball_positions(dims, elements_count, params.equilibrium_distance);

	auto* new_cell = create_cell();
	const index_t cell_index = cells->size() - 1;

	// the cell sits at the centre of mass of its sub-agents
	std::ranges::copy(centre, new_cell->position().begin());
	for (index_t i = 0; i < elements_count; ++i)
	{
		auto position = create_agent(cell_index)->position();
		for (index_t d = 0; d < dims; ++d)
		{
			position[d] = centre[d] + offsets[i * dims + d];
			new_cell->position()[d] += offsets[i * dims + d] / static_cast<real_t>(elements_count);
		}
	}

	return new_cell;
}

namespace {
// Indices of the sub-agents of one cell, in container order
std::vector<index_t> members_of(const agent_data& agents, index_t cell_index)
{
	std::vector<index_t> members;
	for (index_t i = 0; i < agents.agents_count; ++i)
	{
		if (agents.cell_indices[i] == cell_index)
		{
			members.push_back(i);
		}
	}
	return members;
}

// Mean position of the given sub-agents
std::array<real_t, 3> centre_of(const agent_data& agents, const std::vector<index_t>& members)
{
	const index_t dims = agents.base_data.dims;

	std::array<real_t, 3> centre {};
	for (const index_t i : members)
	{
		for (index_t d = 0; d < dims; ++d)
		{
			centre.at(d) += agents.base_data.positions[i * dims + d] / static_cast<real_t>(members.size());
		}
	}
	return centre;
}

// Growth uses its own stream of random numbers, unrelated to the one of the thermal noise
constexpr std::uint64_t growth_stream = 0x67726F7774680000ULL;

// Candidate positions tried at most for a new sub-agent
constexpr index_t growth_candidates = 64;

// A candidate with at least this much room, in equilibrium distances, is taken without looking further
constexpr real_t enough_room = 0.95;
} // namespace

agent_interface* environment::grow_cell(index_t cell_index)
{
	if (cell_index >= cells->size())
	{
		throw std::out_of_range("Invalid cell index");
	}

	const index_t dims = mesh.dims;
	const std::uint64_t seed = mix(params.random_seed ^ growth_stream);
	const index_t growth_step = growth_events++;

	std::array<real_t, 3> best_position {};
	{
		const agent_data& agent_values = generic_agent_solver<agent>::retrieve_agent_data(*agents);
		const std::vector<index_t> members = members_of(agent_values, cell_index);
		const std::vector<real_t>& positions = agent_values.base_data.positions;

		if (members.empty())
		{
			// the first sub-agent of a cell goes where the cell is
			std::ranges::copy(cells->get_agent_at(cell_index)->position(), best_position.begin());
		}

		real_t best_room = -1.0;
		const real_t squared_enough_room = std::pow(enough_room * params.equilibrium_distance, 2);
		for (index_t trial = 0; trial < growth_candidates && !members.empty() && best_room < squared_enough_room;
			 ++trial)
		{
			// one equilibrium distance from a random sub-agent of the cell, in a random direction
			const index_t anchor = members[random_bits(seed, growth_step, trial, dims, 0) % members.size()];

			std::array<real_t, 3> direction {};
			real_t squared_length = 0.0;
			for (index_t d = 0; d < dims; ++d)
			{
				direction.at(d) = standard_normal(seed, growth_step, trial, d);
				squared_length += direction.at(d) * direction.at(d);
			}
			if (squared_length == 0.0)
			{
				direction[0] = 1.0;
				squared_length = 1.0;
			}

			std::array<real_t, 3> position {};
			for (index_t d = 0; d < dims; ++d)
			{
				position.at(d) = positions[anchor * dims + d]
								 + params.equilibrium_distance * direction.at(d) / std::sqrt(squared_length);
			}

			// room: squared distance to the nearest sub-agent of any cell
			real_t room = std::numeric_limits<real_t>::max();
			for (index_t i = 0; i < agent_values.agents_count; ++i)
			{
				real_t squared_distance = 0.0;
				for (index_t d = 0; d < dims; ++d)
				{
					squared_distance +=
						(positions[i * dims + d] - position.at(d)) * (positions[i * dims + d] - position.at(d));
				}
				room = std::min(room, squared_distance);
			}

			if (room > best_room)
			{
				best_room = room;
				best_position = position;
			}
		}
	}

	auto* new_agent = create_agent(cell_index);
	std::ranges::copy_n(best_position.begin(), static_cast<std::ptrdiff_t>(dims), new_agent->position().begin());

	return new_agent;
}

cell_interface* environment::divide_cell(index_t cell_index, std::span<const real_t> normal, real_t daughter_fraction)
{
	if (cell_index >= cells->size())
	{
		throw std::out_of_range("Invalid cell index");
	}

	const index_t dims = mesh.dims;
	if (normal.size() != dims)
	{
		throw std::invalid_argument("The normal of the division plane needs one value per dimension");
	}
	if (std::ranges::all_of(normal, [](real_t x) { return x == 0.0; }))
	{
		throw std::invalid_argument("The normal of the division plane must not be zero");
	}
	// written as a positive condition so that a NaN fraction is rejected as well
	const bool valid_fraction = daughter_fraction > 0.0 && daughter_fraction < 1.0;
	if (!valid_fraction)
	{
		throw std::invalid_argument("The daughter fraction must be between 0 and 1");
	}

	std::vector<index_t> members = members_of(generic_agent_solver<agent>::retrieve_agent_data(*agents), cell_index);
	if (members.size() < 2)
	{
		throw std::invalid_argument("A cell needs at least two sub-agents to divide");
	}

	// order the sub-agents along the normal; equal projections keep their container order
	{
		const std::vector<real_t>& positions =
			generic_agent_solver<agent>::retrieve_agent_data(*agents).base_data.positions;
		const auto projection = [&](index_t i) {
			real_t sum = 0.0;
			for (index_t d = 0; d < dims; ++d)
			{
				sum += positions[i * dims + d] * normal[d];
			}
			return sum;
		};
		std::ranges::stable_sort(members, [&](index_t a, index_t b) { return projection(a) < projection(b); });
	}

	// the daughter takes the sub-agents furthest along the normal; both cells keep at least one
	const auto count = static_cast<real_t>(members.size());
	const auto daughter_count =
		std::clamp<index_t>(static_cast<index_t>(std::llround(daughter_fraction * count)), 1, members.size() - 1);
	const index_t parent_count = members.size() - daughter_count;

	auto* daughter = create_cell();
	const index_t daughter_index = cells->size() - 1;
	auto* parent = cells->get_agent_at(cell_index);

	daughter->is_movable() = parent->is_movable();
	daughter->is_motile() = parent->is_motile();
	daughter->migration_speed() = parent->migration_speed();
	daughter->migration_bias() = parent->migration_bias();
	daughter->persistence_time() = parent->persistence_time();
	std::ranges::copy(parent->migration_bias_direction(), daughter->migration_bias_direction().begin());
	std::ranges::copy(parent->motility_direction(), daughter->motility_direction().begin());
	std::ranges::copy(parent->velocity(), daughter->velocity().begin());

	agent_data& agent_values = generic_agent_solver<agent>::retrieve_agent_data(*agents);
	cell_data& cell_values = generic_agent_solver<cell>::retrieve_agent_data(*cells);

	const std::vector<index_t> staying(members.begin(), members.begin() + static_cast<std::ptrdiff_t>(parent_count));
	const std::vector<index_t> leaving(members.begin() + static_cast<std::ptrdiff_t>(parent_count), members.end());

	for (const index_t i : leaving)
	{
		agent_values.cell_indices[i] = daughter_index;
	}
	cell_values.agent_counts[cell_index] = parent_count;
	cell_values.agent_counts[daughter_index] = daughter_count;

	// both cells sit at the centres of mass of their sub-agents
	const std::array<real_t, 3> parent_centre = centre_of(agent_values, staying);
	const std::array<real_t, 3> daughter_centre = centre_of(agent_values, leaving);
	std::ranges::copy_n(parent_centre.begin(), static_cast<std::ptrdiff_t>(dims), parent->position().begin());
	std::ranges::copy_n(daughter_centre.begin(), static_cast<std::ptrdiff_t>(dims), daughter->position().begin());

	return daughter;
}

std::vector<real_t> environment::longest_axis(index_t cell_index)
{
	if (cell_index >= cells->size())
	{
		throw std::out_of_range("Invalid cell index");
	}

	const index_t dims = mesh.dims;
	const agent_data& agent_values = generic_agent_solver<agent>::retrieve_agent_data(*agents);
	const std::vector<index_t> members = members_of(agent_values, cell_index);

	std::vector<real_t> axis(dims, 0.0);
	axis[0] = 1.0;
	if (members.size() < 2)
	{
		return axis;
	}

	// covariance of the sub-agent positions
	const std::array<real_t, 3> centre = centre_of(agent_values, members);
	std::array<std::array<real_t, 3>, 3> covariance {};
	for (const index_t i : members)
	{
		for (index_t a = 0; a < dims; ++a)
		{
			for (index_t b = 0; b < dims; ++b)
			{
				covariance.at(a).at(b) += (agent_values.base_data.positions[i * dims + a] - centre.at(a))
										  * (agent_values.base_data.positions[i * dims + b] - centre.at(b));
			}
		}
	}

	// Its eigenvector with the largest eigenvalue, by power iteration. The start is the coordinate
	// axis with the largest spread, which cannot be perpendicular to that eigenvector.
	index_t widest = 0;
	for (index_t d = 1; d < dims; ++d)
	{
		if (covariance.at(d).at(d) > covariance.at(widest).at(widest))
		{
			widest = d;
		}
	}
	if (covariance.at(widest).at(widest) == 0.0)
	{
		return axis;
	}
	std::ranges::fill(axis, 0.0);
	axis[widest] = 1.0;

	for (int iteration = 0; iteration < 200; ++iteration)
	{
		std::array<real_t, 3> next {};
		real_t squared_length = 0.0;
		for (index_t a = 0; a < dims; ++a)
		{
			for (index_t b = 0; b < dims; ++b)
			{
				next.at(a) += covariance.at(a).at(b) * axis[b];
			}
			squared_length += next.at(a) * next.at(a);
		}
		for (index_t d = 0; d < dims; ++d)
		{
			axis[d] = next.at(d) / std::sqrt(squared_length);
		}
	}

	// an axis has no sign: by convention its largest component is positive
	const auto largest = std::ranges::max_element(axis, {}, [](real_t x) { return std::abs(x); });
	if (*largest < 0.0)
	{
		for (real_t& x : axis)
		{
			x = -x;
		}
	}

	return axis;
}

agent_interface* environment::create_agent(index_t cell_index)
{
	if (cell_index >= cells->size())
	{
		throw std::out_of_range("Invalid cell index");
	}

	auto* new_agent = agents->create();
	new_agent->cell_index() = cell_index;

	++generic_agent_solver<cell>::retrieve_agent_data(*cells).agent_counts[cell_index];

	return new_agent;
}

void environment::remove_agent(index_t agent_index)
{
	if (agent_index >= agents->size())
	{
		throw std::out_of_range("Invalid agent index");
	}

	const index_t cell_index = generic_agent_solver<agent>::retrieve_agent_data(*agents).cell_indices[agent_index];
	if (cell_index != no_cell)
	{
		--generic_agent_solver<cell>::retrieve_agent_data(*cells).agent_counts[cell_index];
	}

	agents->remove_at(agent_index);
}

void environment::remove_cell(index_t cell_index)
{
	if (cell_index >= cells->size())
	{
		throw std::out_of_range("Invalid cell index");
	}

	auto& cell_indices = generic_agent_solver<agent>::retrieve_agent_data(*agents).cell_indices;

	// Going backwards, the sub-agent moved into a freed slot has already been checked.
	for (index_t i = agents->size(); i-- > 0;)
	{
		if (cell_indices[i] == cell_index)
		{
			agents->remove_at(i);
		}
	}

	// The last cell is about to take over the freed index.
	const index_t last_cell_index = cells->size() - 1;
	if (cell_index != last_cell_index)
	{
		std::ranges::replace(cell_indices, last_cell_index, cell_index);
	}

	cells->remove_at(cell_index);
}

} // namespace physicore::mechanics::micromechanics
