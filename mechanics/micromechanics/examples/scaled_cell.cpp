// Builds one cell from cell-level parameters with a chosen number of elements, lets it relax and prints
// how its size develops, as CSV.
//
// The cell is the one of Sandersius and Newman (2008), section 4.1, in micrometres, seconds and
// piconewtons: radius 10, stiffness 5000, damping 5000. The header line lists the element parameters
// that scale_to_elements() derived, so documentation can be checked against them.
//
// Usage: mechanics.micromechanics.scaled_cell [elements_per_cell [steepness [end_time [snapshots]]]]

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <span>
#include <string>

#include <micromechanics/cell_parameters.h>
#include <micromechanics/environment.h>
#include <micromechanics/solver_registry.h>

using namespace physicore;
using namespace physicore::mechanics::micromechanics;

namespace {

struct cell_size
{
	// sqrt(5/3) times the radius of gyration: the radius of a uniform ball with the same spread
	real_t radius_from_gyration = 0.0;
	// largest distance of an element centre from the centre of the cell
	real_t largest_distance = 0.0;
	real_t largest_speed = 0.0;
};

cell_size measure(environment& env)
{
	const index_t count = env.agents->size();

	std::array<real_t, 3> centre {};
	for (index_t i = 0; i < count; ++i)
	{
		const auto position = env.agents->get_agent_at(i)->position();
		for (std::size_t d = 0; d < 3; ++d)
		{
			centre.at(d) += position[d] / static_cast<real_t>(count);
		}
	}

	cell_size size;
	real_t squared_gyration_radius = 0.0;
	for (index_t i = 0; i < count; ++i)
	{
		auto* element = env.agents->get_agent_at(i);
		const auto position = element->position();

		real_t squared_distance = 0.0;
		for (std::size_t d = 0; d < 3; ++d)
		{
			squared_distance += std::pow(position[d] - centre.at(d), 2);
			size.largest_speed = std::max(size.largest_speed, std::abs(element->velocity()[d]));
		}
		squared_gyration_radius += squared_distance / static_cast<real_t>(count);
		size.largest_distance = std::max(size.largest_distance, std::sqrt(squared_distance));
	}
	size.radius_from_gyration = std::sqrt(5.0 / 3.0 * squared_gyration_radius);

	return size;
}

} // namespace

int main(int argc, char** argv)
{
	const std::span<char*> args(argv, static_cast<std::size_t>(argc));

	index_t elements_per_cell = 256;
	real_t end_time = 6.0;
	index_t snapshots = 24;
	simulation_parameters shape;

	try
	{
		if (args.size() > 5)
		{
			throw std::invalid_argument("expected at most 4 arguments");
		}
		if (args.size() > 1)
		{
			elements_per_cell = std::stoull(args[1]);
		}
		if (args.size() > 2)
		{
			shape.steepness = std::stod(args[2]);
		}
		if (args.size() > 3)
		{
			end_time = std::stod(args[3]);
		}
		if (args.size() > 4)
		{
			snapshots = std::stoull(args[4]);
		}
		if (end_time <= 0.0 || snapshots == 0)
		{
			throw std::invalid_argument("end_time and snapshots must be positive");
		}

		const cell_parameters cell { .radius = 10.0, .stiffness = 5000.0, .damping = 5000.0 };
		const cell_type type = scale_to_elements(cell, elements_per_cell, shape);

		// a twentieth of the largest stable timestep of an isolated pair
		const real_t pair_stiffness = 8.0 * type.intra_well_depth * shape.steepness * shape.steepness
									  / (type.equilibrium_distance * type.equilibrium_distance);
		const real_t timestep = 0.05 * type.drag_coefficient / pair_stiffness;

		environment env({ 3, { -100, -100, -100 }, { 100, 100, 100 }, { 20, 20, 20 } }, timestep);
		env.params = shape;
		env.cell_types[0] = type;
		env.create_cell();

		// elements scattered uniformly in a ball of the requested radius
		std::mt19937 generator(1);
		std::uniform_real_distribution<real_t> coordinate(-1.0, 1.0);
		while (env.agents->size() < elements_per_cell)
		{
			const real_t x = coordinate(generator);
			const real_t y = coordinate(generator);
			const real_t z = coordinate(generator);
			if (x * x + y * y + z * z <= 1.0)
			{
				auto position = env.create_agent(0)->position();
				position[0] = cell.radius * x;
				position[1] = cell.radius * y;
				position[2] = cell.radius * z;
			}
		}

		env.solver = solver_registry::instance().get("openmp_solver");
		if (!env.solver)
		{
			throw std::runtime_error("the openmp_solver is not registered");
		}
		env.solver->initialize(env);

		const auto steps_per_snapshot = std::max<index_t>(
			1, static_cast<index_t>(std::llround(end_time / (timestep * static_cast<real_t>(snapshots)))));

		std::cout << std::setprecision(std::numeric_limits<real_t>::max_digits10);
		std::cout << "# elements_per_cell=" << elements_per_cell << " cell_radius=" << cell.radius
				  << " cell_stiffness=" << cell.stiffness << " cell_damping=" << cell.damping
				  << " packing_density=" << cell.packing_density << " scaling_correction=" << cell.scaling_correction
				  << " steepness=" << shape.steepness << " equilibrium_distance=" << type.equilibrium_distance
				  << " drag_coefficient=" << type.drag_coefficient << " intra_well_depth=" << type.intra_well_depth
				  << " pair_stiffness=" << pair_stiffness << " timestep=" << timestep << '\n';
		std::cout << "time,radius_from_gyration,largest_distance,largest_speed\n";

		for (index_t snapshot = 0; snapshot <= snapshots; ++snapshot)
		{
			if (snapshot > 0)
			{
				env.solver->solve(env, steps_per_snapshot);
			}
			const cell_size size = measure(env);
			std::cout << static_cast<real_t>(snapshot * steps_per_snapshot) * timestep << ','
					  << size.radius_from_gyration << ',' << size.largest_distance << ',' << size.largest_speed << '\n';
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << "scaled_cell: " << e.what() << '\n'
				  << "usage: " << args[0] << " [elements_per_cell [steepness [end_time [snapshots]]]]\n";
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
