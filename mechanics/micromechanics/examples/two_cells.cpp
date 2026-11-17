// Two cells of 27 sub-agents each are placed side by side and left to relax. Prints the positions of
// all sub-agents at regular times as CSV.
//
// The outcome depends on the ratio of the inter- to the intra-cellular well depth: the cells ignore
// each other (0), adhere without mixing (around 0.5) or merge into one aggregate (1).
//
// Usage: mechanics.micromechanics.two_cells [inter_to_intra_ratio [end_time [snapshots [vtk_output_dir]]]]
//
// With vtk_output_dir, every snapshot is also written as a VTK file; open
// <vtk_output_dir>/micromechanics_agents.pvd in ParaView.

#include <cmath>
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <span>
#include <string>

#include <micromechanics/environment.h>
#include <micromechanics/solver_registry.h>
#include <micromechanics/vtk_serializer.h>

using namespace physicore;
using namespace physicore::mechanics::micromechanics;

namespace {

// One cell on a 3 x 3 x 3 lattice with unit spacing, every sub-agent displaced by up to +-0.15 per axis
void add_cell(environment& env, real_t x_offset, unsigned seed)
{
	std::mt19937 generator(seed);
	std::uniform_real_distribution<real_t> displacement(-0.15, 0.15);

	env.create_cell();
	const index_t cell_index = env.cells->size() - 1;

	for (int x = 0; x < 3; ++x)
	{
		for (int y = 0; y < 3; ++y)
		{
			for (int z = 0; z < 3; ++z)
			{
				auto position = env.create_agent(cell_index)->position();
				position[0] = x_offset + x + displacement(generator);
				position[1] = y + displacement(generator);
				position[2] = z + displacement(generator);
			}
		}
	}
}

void print_snapshot(environment& env, real_t time)
{
	for (index_t i = 0; i < env.agents->size(); ++i)
	{
		auto* element = env.agents->get_agent_at(i);
		const auto position = element->position();

		std::cout << time << ',' << i << ',' << element->cell_index() << ',' << position[0] << ',' << position[1] << ','
				  << position[2] << '\n';
	}
}

} // namespace

int main(int argc, char** argv)
{
	const std::span<char*> args(argv, static_cast<std::size_t>(argc));

	real_t ratio = 0.5;
	real_t end_time = 10.0;
	index_t snapshots = 20;
	std::string vtk_output_dir;

	try
	{
		if (args.size() > 5)
		{
			throw std::invalid_argument("expected at most 4 arguments");
		}
		if (args.size() > 1)
		{
			ratio = std::stod(args[1]);
		}
		if (args.size() > 2)
		{
			end_time = std::stod(args[2]);
		}
		if (args.size() > 3)
		{
			snapshots = std::stoull(args[3]);
		}
		if (args.size() > 4)
		{
			vtk_output_dir = args[4];
		}
		if (ratio < 0.0 || end_time <= 0.0 || snapshots == 0)
		{
			throw std::invalid_argument("ratio must not be negative; end_time and snapshots must be positive");
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << "two_cells: " << e.what() << '\n'
				  << "usage: " << args[0] << " [inter_to_intra_ratio [end_time [snapshots [vtk_output_dir]]]]\n";
		return EXIT_FAILURE;
	}

	const real_t timestep = 0.002;
	environment env({ 3, { -100, -100, -100 }, { 100, 100, 100 }, { 20, 20, 20 } }, timestep);
	env.params.inter_well_depth = ratio * env.params.intra_well_depth;

	// the second cell starts 1.3 equilibrium distances away from the first, inside the interaction range
	add_cell(env, 0.0, 42);
	add_cell(env, 3.3, 43);

	env.solver = solver_registry::instance().get("openmp_solver");
	if (!env.solver)
	{
		std::cerr << "two_cells: the openmp_solver is not registered\n";
		return EXIT_FAILURE;
	}
	env.solver->initialize(env);

	if (!vtk_output_dir.empty())
	{
		env.serializer = make_vtk_agents_serializer(vtk_output_dir);
	}

	const auto steps_per_snapshot =
		static_cast<index_t>(std::llround(end_time / (timestep * static_cast<real_t>(snapshots))));

	std::cout << std::setprecision(std::numeric_limits<real_t>::max_digits10);
	std::cout << "# inter_to_intra_ratio=" << ratio << " timestep=" << timestep
			  << " equilibrium_distance=" << env.params.equilibrium_distance << '\n';
	std::cout << "time,agent,cell,x,y,z\n";

	print_snapshot(env, 0.0);
	env.serialize_state(0.0);
	for (index_t snapshot = 1; snapshot <= snapshots; ++snapshot)
	{
		env.solver->solve(env, steps_per_snapshot);
		const real_t time = static_cast<real_t>(snapshot * steps_per_snapshot) * timestep;
		print_snapshot(env, time);
		env.serialize_state(time);
	}

	return EXIT_SUCCESS;
}
