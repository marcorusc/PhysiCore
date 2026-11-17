// A single cell grows and divides for a few generations. In every generation each cell doubles its
// number of sub-agents, one at a time, and then divides across its longest axis.
//
// Prints the positions of all sub-agents after every growth phase and after every division as CSV.
// With daughter_fraction 0.5 all cells stay the same size; any other value gives a population of
// cells of different sizes.
//
// Usage: mechanics.micromechanics.growth_division [daughter_fraction [generations [birth_count
//                                                  [steepness [vtk_output_dir]]]]]

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <span>
#include <string>
#include <vector>

#include <micromechanics/environment.h>
#include <micromechanics/solver_registry.h>
#include <micromechanics/vtk_serializer.h>

using namespace physicore;
using namespace physicore::mechanics::micromechanics;

namespace {

struct options
{
	real_t daughter_fraction = 0.5;
	index_t generations = 3;
	index_t birth_count = 32;
	real_t steepness = 2.0;
	std::string vtk_output_dir;
};

options parse(std::span<char*> args)
{
	if (args.size() > 6)
	{
		throw std::invalid_argument("expected at most 5 arguments");
	}

	options result;
	if (args.size() > 1)
	{
		result.daughter_fraction = std::stod(args[1]);
	}
	if (args.size() > 2)
	{
		result.generations = std::stoull(args[2]);
	}
	if (args.size() > 3)
	{
		result.birth_count = std::stoull(args[3]);
	}
	if (args.size() > 4)
	{
		result.steepness = std::stod(args[4]);
	}
	if (args.size() > 5)
	{
		result.vtk_output_dir = args[5];
	}
	if (result.birth_count < 2 || result.steepness <= 0.0)
	{
		throw std::invalid_argument("birth_count must be at least 2 and steepness positive");
	}
	return result;
}

// The simulation together with its clock
struct simulation
{
	environment env;
	real_t timestep;
	real_t time = 0.0;

	index_t steps(real_t duration) const { return static_cast<index_t>(std::llround(duration / timestep)); }

	// Advances by the given number of quarter time units, writing a VTK frame after each
	void relax(index_t quarters)
	{
		for (index_t quarter = 0; quarter < quarters; ++quarter)
		{
			env.solver->solve(env, steps(0.25));
			time += 0.25;
			env.serialize_state(time);
		}
	}

	void print_snapshot(const std::string& phase, index_t generation)
	{
		for (index_t i = 0; i < env.agents->size(); ++i)
		{
			auto* element = env.agents->get_agent_at(i);
			const auto position = element->position();

			std::cout << phase << ',' << generation << ',' << time << ',' << i << ',' << element->cell_index() << ','
					  << position[0] << ',' << position[1] << ',' << position[2] << '\n';
		}
	}

	// Every cell doubles the number of sub-agents it has now. One sub-agent per cell in turn, with a
	// short relaxation after each round, so that all cells grow at the same time.
	void grow_all_cells()
	{
		const index_t cells_count = env.cells->size();

		std::vector<index_t> to_add(cells_count);
		for (index_t c = 0; c < cells_count; ++c)
		{
			to_add[c] = env.cells->get_agent_at(c)->agent_count();
		}
		const index_t rounds = *std::ranges::max_element(to_add);

		for (index_t round = 0; round < rounds; ++round)
		{
			for (index_t c = 0; c < cells_count; ++c)
			{
				if (round < to_add[c])
				{
					env.grow_cell(c);
				}
			}
			env.solver->solve(env, steps(0.05));
			time += 0.05;
		}
	}

	// Every cell divides across its longest axis
	void divide_all_cells(real_t daughter_fraction)
	{
		const index_t cells_count = env.cells->size();
		for (index_t c = 0; c < cells_count; ++c)
		{
			if (env.cells->get_agent_at(c)->agent_count() >= 2)
			{
				env.divide_cell(c, env.longest_axis(c), daughter_fraction);
			}
		}
	}
};

void run(const options& opts)
{
	simulation_parameters params;
	params.steepness = opts.steepness;
	const cell_type type;

	// a twentieth of the largest stable timestep of an isolated pair
	const real_t pair_stiffness = 8.0 * type.intra_well_depth * params.steepness * params.steepness
								  / (type.equilibrium_distance * type.equilibrium_distance);
	const real_t timestep = 0.05 * type.drag_coefficient / pair_stiffness;

	simulation sim { .env = { { 3, { -100, -100, -100 }, { 100, 100, 100 }, { 20, 20, 20 } }, timestep },
					 .timestep = timestep };
	sim.env.params = params;

	sim.env.solver = solver_registry::instance().get("openmp_solver");
	if (!sim.env.solver)
	{
		throw std::runtime_error("the openmp_solver is not registered");
	}
	if (!opts.vtk_output_dir.empty())
	{
		sim.env.serializer = make_vtk_agents_serializer(opts.vtk_output_dir);
	}

	const std::array<real_t, 3> centre = { 0.0, 0.0, 0.0 };
	sim.env.create_cell(centre, opts.birth_count);
	sim.env.solver->initialize(sim.env);

	std::cout << std::setprecision(8);
	std::cout << "# daughter_fraction=" << opts.daughter_fraction << " generations=" << opts.generations
			  << " birth_count=" << opts.birth_count << " steepness=" << params.steepness << " timestep=" << timestep
			  << '\n';
	std::cout << "phase,generation,time,agent,cell,x,y,z\n";

	sim.env.serialize_state(sim.time);
	sim.relax(4);
	sim.print_snapshot("start", 0);

	for (index_t generation = 1; generation <= opts.generations; ++generation)
	{
		sim.grow_all_cells();
		sim.relax(4);
		sim.print_snapshot("grown", generation);

		sim.divide_all_cells(opts.daughter_fraction);
		sim.relax(8);
		sim.print_snapshot("divided", generation);
	}
}

} // namespace

int main(int argc, char** argv)
{
	const std::span<char*> args(argv, static_cast<std::size_t>(argc));

	try
	{
		run(parse(args));
	}
	catch (const std::exception& e)
	{
		std::cerr << "growth_division: " << e.what() << '\n'
				  << "usage: " << args[0]
				  << " [daughter_fraction [generations [birth_count [steepness [vtk_output_dir]]]]]\n";
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
