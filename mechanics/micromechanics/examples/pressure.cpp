// Two small experiments on the pressure of a cell, printed as CSV.
//
//   pair [elements_per_cell]
//       Two cells in three dimensions, held at a series of centre distances. Nothing moves. Prints the
//       pressure on one cell against the distance: positive where the cells overlap, negative where
//       they only attract, zero beyond the cutoff.
//
//   confined [spacing]
//       Nineteen cells in two dimensions on a hexagonal layout: one in the middle, a ring of six and
//       a ring of twelve. The outer ring cannot move and holds the other seven in place, which then
//       relax. Prints the pressure of every cell at the end. A smaller spacing squeezes the cells.
//
// Usage: mechanics.micromechanics.pressure pair|confined [value]

#include <array>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <span>
#include <stdexcept>
#include <string>

#include <micromechanics/environment.h>
#include <micromechanics/solver_registry.h>

using namespace physicore;
using namespace physicore::mechanics::micromechanics;

namespace {

environment make_environment(index_t dims, real_t timestep)
{
	return { { dims, { -100, -100, -100 }, { 100, 100, 100 }, { 20, 20, 20 } }, timestep };
}

void attach_solver(environment& env)
{
	env.solver = solver_registry::instance().get("openmp_solver");
	if (!env.solver)
	{
		throw std::runtime_error("the openmp_solver is not registered");
	}
	env.solver->initialize(env);
}

void run_pair(index_t elements_per_cell)
{
	const simulation_parameters params;

	std::cout << "# experiment=pair elements_per_cell=" << elements_per_cell
			  << " equilibrium_distance=" << params.equilibrium_distance
			  << " inter_well_depth=" << params.inter_well_depth << " cutoff_factor=" << params.cutoff_factor << '\n';
	std::cout << "distance,pressure\n";

	for (int step = 0; step <= 160; ++step)
	{
		const real_t distance = 1.0 + 0.05 * step;

		environment env = make_environment(3, 0.001);
		for (const real_t side : { -0.5, 0.5 })
		{
			const std::array<real_t, 3> centre = { side * distance, 0.0, 0.0 };
			env.create_cell(centre, elements_per_cell);
		}
		// initialising the solver computes the pressure of the cells as they were placed
		attach_solver(env);

		std::cout << distance << ',' << env.cells->get_agent_at(0)->pressure() << '\n';
	}
}

void run_confined(real_t spacing)
{
	const index_t elements_per_cell = 37; // a hexagonal disc that reaches 3 from its centre
	const real_t timestep = 0.002;
	const index_t steps = 10000;

	environment env = make_environment(2, timestep);

	// hexagonal layout in axial coordinates: all sites within two steps of the middle
	for (int q = -2; q <= 2; ++q)
	{
		for (int r = -2; r <= 2; ++r)
		{
			const int ring = (std::abs(q) + std::abs(r) + std::abs(q + r)) / 2;
			if (ring > 2)
			{
				continue;
			}
			const std::array<real_t, 2> centre = { spacing * (q + 0.5 * r), spacing * 0.5 * std::numbers::sqrt3 * r };
			auto* cell = env.create_cell(centre, elements_per_cell);
			cell->is_movable() = ring == 2 ? 0 : 1;
		}
	}
	attach_solver(env);
	env.solver->solve(env, steps);

	std::cout << "# experiment=confined spacing=" << spacing << " elements_per_cell=" << elements_per_cell
			  << " cells=" << env.cells->size() << " time=" << timestep * static_cast<real_t>(steps) << '\n';
	std::cout << "cell,ring,x,y,speed,pressure\n";

	index_t index = 0;
	for (int q = -2; q <= 2; ++q)
	{
		for (int r = -2; r <= 2; ++r)
		{
			const int ring = (std::abs(q) + std::abs(r) + std::abs(q + r)) / 2;
			if (ring > 2)
			{
				continue;
			}
			auto* cell = env.cells->get_agent_at(index);
			std::cout << index << ',' << ring << ',' << cell->position()[0] << ',' << cell->position()[1] << ','
					  << std::hypot(cell->velocity()[0], cell->velocity()[1]) << ',' << cell->pressure() << '\n';
			++index;
		}
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
		if (experiment == "pair")
		{
			const index_t elements = args.size() > 2 ? std::stoull(args[2]) : 55;
			if (elements < 2)
			{
				throw std::invalid_argument("a cell needs at least 2 elements to have a pressure");
			}
			run_pair(elements);
		}
		else if (experiment == "confined")
		{
			const real_t spacing = args.size() > 2 ? std::stod(args[2]) : 6.5;
			if (spacing <= 0.0)
			{
				throw std::invalid_argument("the spacing must be positive");
			}
			run_confined(spacing);
		}
		else
		{
			throw std::invalid_argument("unknown experiment '" + experiment + "'");
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << "pressure: " << e.what() << '\n' << "usage: " << args[0] << " pair|confined [value]\n";
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
