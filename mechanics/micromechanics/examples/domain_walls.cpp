// Two small experiments on the walls of the domain, printed as CSV.
//
//   force
//       The force and the energy of a wall against the distance of an element from it, next to the
//       force between two elements at twice that distance.
//
//   stall [migration_speed]
//       A motile cell of 37 elements in two dimensions crawls into a wall. Prints over time how far
//       its centre is from the wall, its speed, its pressure and how often the hard stop of the
//       solver was needed.
//
// Usage: mechanics.micromechanics.domain_walls force|stall [migration_speed]

#include <array>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <span>
#include <stdexcept>
#include <string>

#include <micromechanics/domain_wall.h>
#include <micromechanics/environment.h>
#include <micromechanics/solver_registry.h>

using namespace physicore;
using namespace physicore::mechanics::micromechanics;

namespace {

void run_force()
{
	const simulation_parameters params;
	const pair_potential wall = params.wall_potential();

	std::cout << "# experiment=force wall_well_depth=" << params.wall_well_depth
			  << " equilibrium_distance=" << params.equilibrium_distance << " steepness=" << params.steepness << '\n';
	std::cout << "distance,wall_force,wall_energy,pair_force_at_twice_the_distance\n";

	for (int step = 0; step <= 200; ++step)
	{
		const real_t distance = 0.005 * step;
		// a wall at 0, the other one far away
		std::cout << distance << ',' << wall_force(wall, distance, 0.0, 1000.0) << ','
				  << wall_energy(wall, distance, 0.0, 1000.0) << ','
				  << force_coefficient(wall, 4.0 * distance * distance) * 2.0 * distance << '\n';
	}
}

void run_stall(real_t speed)
{
	const index_t elements_per_cell = 37;
	const real_t timestep = 0.002;
	const index_t steps_per_snapshot = 50;
	const real_t wall_position = 12.0;

	environment env({ 2, { -12, -12, -12 }, { 12, 12, 12 }, { 1, 1, 1 } }, timestep);
	env.params.domain_walls = true;

	const std::array<real_t, 2> centre = { 4.0, 0.0 };
	auto* cell = env.create_cell(centre, elements_per_cell);
	cell->is_motile() = 1;
	cell->migration_speed() = speed;
	cell->persistence_time() = 1.0;
	cell->migration_bias() = 1.0;
	cell->migration_bias_direction()[0] = 1.0; // towards the wall at x = 12

	env.solver = solver_registry::instance().get("openmp_solver");
	if (!env.solver)
	{
		throw std::runtime_error("the openmp_solver is not registered");
	}
	env.solver->initialize(env);

	const real_t end_time = 8.0 / speed + 15.0; // time to reach the wall, and then to settle
	const auto snapshots = static_cast<index_t>(std::llround(end_time / (timestep * steps_per_snapshot)));

	std::cout << "# experiment=stall migration_speed=" << speed << " elements_per_cell=" << elements_per_cell
			  << " drag_coefficient=" << env.params.drag_coefficient
			  << " equilibrium_distance=" << env.params.equilibrium_distance << " timestep=" << timestep << '\n';
	std::cout << "time,distance_to_wall,nearest_element_to_wall,speed,pressure\n";

	for (index_t snapshot = 0; snapshot <= snapshots; ++snapshot)
	{
		if (snapshot > 0)
		{
			env.solver->solve(env, steps_per_snapshot);
		}

		real_t nearest = wall_position;
		for (index_t i = 0; i < env.agents->size(); ++i)
		{
			nearest = std::min(nearest, wall_position - env.agents->get_agent_at(i)->position()[0]);
		}

		std::cout << static_cast<real_t>(snapshot * steps_per_snapshot) * timestep << ','
				  << wall_position - cell->position()[0] << ',' << nearest << ','
				  << std::hypot(cell->velocity()[0], cell->velocity()[1]) << ',' << cell->pressure() << '\n';
	}
}

} // namespace

int main(int argc, char** argv)
{
	const std::span<char*> args(argv, static_cast<std::size_t>(argc));

	try
	{
		if (args.size() < 2 || args.size() > 3)
		{
			throw std::invalid_argument("expected an experiment and at most one value");
		}

		std::cout << std::setprecision(9);

		const std::string experiment = args[1];
		if (experiment == "force" && args.size() == 2)
		{
			run_force();
		}
		else if (experiment == "stall")
		{
			const real_t speed = args.size() > 2 ? std::stod(args[2]) : 0.6;
			if (speed <= 0.0)
			{
				throw std::invalid_argument("the migration speed must be positive");
			}
			run_stall(speed);
		}
		else
		{
			throw std::invalid_argument("unknown experiment '" + experiment + "'");
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << "domain_walls: " << e.what() << '\n' << "usage: " << args[0] << " force|stall [migration_speed]\n";
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
