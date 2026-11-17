// Two small experiments on cell motility (a persistent random walk), printed as CSV.
//
//   msd [dims [persistence_time [bias]]]
//       Many motile cells that do not interact. Prints over time the mean squared displacement and
//       the mean displacement along the bias direction (the x axis). Without bias the first should
//       follow 2 v^2 tau (t - tau (1 - exp(-t / tau))) with v the migration speed and tau the
//       persistence time.
//
//   tracks [dims [persistence_time [bias]]]
//       A few motile cells of 13 sub-agents each, far apart. Prints the position of every cell over
//       time.
//
// The migration speed is 1 and the bias direction the x axis.
//
// Usage: mechanics.micromechanics.motility msd|tracks [dims [persistence_time [bias]]]

#include <array>
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include <micromechanics/environment.h>
#include <micromechanics/solver_registry.h>

using namespace physicore;
using namespace physicore::mechanics::micromechanics;

namespace {

constexpr real_t migration_speed = 1.0;
constexpr real_t timestep = 0.01;
constexpr real_t separation = 20.0;

struct walk_parameters
{
	index_t dims = 2;
	real_t persistence_time = 1.0;
	real_t bias = 0.0;
};

environment make_environment(index_t dims)
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

void make_motile(cell_interface& cell, const walk_parameters& walk)
{
	cell.is_motile() = 1;
	cell.migration_speed() = migration_speed;
	cell.persistence_time() = walk.persistence_time;
	cell.migration_bias() = walk.bias;
	cell.migration_bias_direction()[0] = 1.0;
}

void print_header(const std::string& experiment, const walk_parameters& walk, index_t cells)
{
	std::cout << "# experiment=" << experiment << " dims=" << walk.dims << " migration_speed=" << migration_speed
			  << " persistence_time=" << walk.persistence_time << " bias=" << walk.bias << " cells=" << cells
			  << " timestep=" << timestep << '\n';
}

void run_msd(const walk_parameters& walk)
{
	const index_t cells = 2000;
	const index_t snapshots = 80;
	const index_t steps_per_snapshot = 10;

	environment env = make_environment(walk.dims);
	env.params.random_seed = 2026;
	// one sub-agent per cell and no force between cells: the cells are free
	env.params.inter_well_depth = 0.0;
	for (index_t c = 0; c < cells; ++c)
	{
		make_motile(*env.create_cell(), walk);
		env.create_agent(c)->position()[0] = separation * static_cast<real_t>(c);
	}
	attach_solver(env);

	print_header("msd", walk, cells);
	std::cout << "time,mean_squared_displacement,mean_displacement_along_bias\n";

	for (index_t snapshot = 0; snapshot <= snapshots; ++snapshot)
	{
		if (snapshot > 0)
		{
			env.solver->solve(env, steps_per_snapshot);
		}

		real_t squared = 0.0;
		real_t along_bias = 0.0;
		for (index_t c = 0; c < cells; ++c)
		{
			const auto position = env.cells->get_agent_at(c)->position();
			for (index_t d = 0; d < walk.dims; ++d)
			{
				const real_t start = d == 0 ? separation * static_cast<real_t>(c) : 0.0;
				squared += (position[d] - start) * (position[d] - start);
			}
			along_bias += position[0] - separation * static_cast<real_t>(c);
		}
		std::cout << static_cast<real_t>(snapshot * steps_per_snapshot) * timestep << ','
				  << squared / static_cast<real_t>(cells) << ',' << along_bias / static_cast<real_t>(cells) << '\n';
	}
}

void run_tracks(const walk_parameters& walk)
{
	const index_t cells = 8;
	const index_t elements_per_cell = 13;
	const index_t snapshots = 400;
	const index_t steps_per_snapshot = 5;

	environment env = make_environment(walk.dims);
	env.params.random_seed = 2026;
	// the cells start far apart, but their walks may cross: no force between cells keeps them free
	env.params.inter_well_depth = 0.0;
	for (index_t c = 0; c < cells; ++c)
	{
		const std::array<real_t, 3> centre = { separation * static_cast<real_t>(c), 0.0, 0.0 };
		make_motile(*env.create_cell(std::span<const real_t>(centre.data(), walk.dims), elements_per_cell), walk);
	}
	attach_solver(env);

	print_header("tracks", walk, cells);
	std::cout << "time,cell,x,y,z\n";

	for (index_t snapshot = 0; snapshot <= snapshots; ++snapshot)
	{
		if (snapshot > 0)
		{
			env.solver->solve(env, steps_per_snapshot);
		}

		for (index_t c = 0; c < cells; ++c)
		{
			const auto position = env.cells->get_agent_at(c)->position();
			// every cell is printed relative to its own start
			std::cout << static_cast<real_t>(snapshot * steps_per_snapshot) * timestep << ',' << c << ','
					  << position[0] - separation * static_cast<real_t>(c);
			for (index_t d = 1; d < 3; ++d)
			{
				std::cout << ',' << (d < walk.dims ? position[d] : 0.0);
			}
			std::cout << '\n';
		}
	}
}

} // namespace

int main(int argc, char** argv)
{
	const std::span<char*> args(argv, static_cast<std::size_t>(argc));

	try
	{
		if (args.size() < 2 || args.size() > 5)
		{
			throw std::invalid_argument("expected an experiment and at most 3 further arguments");
		}

		walk_parameters walk;
		if (args.size() > 2)
		{
			walk.dims = std::stoull(args[2]);
		}
		if (args.size() > 3)
		{
			walk.persistence_time = std::stod(args[3]);
		}
		if (args.size() > 4)
		{
			walk.bias = std::stod(args[4]);
		}
		if (walk.dims < 1 || walk.dims > 3 || walk.persistence_time <= 0.0 || walk.bias < 0.0 || walk.bias > 1.0)
		{
			throw std::invalid_argument(
				"dims must be 1, 2 or 3, the persistence time positive and the bias between 0 and 1");
		}

		// nine digits are plenty for a plot and keep the output small
		std::cout << std::setprecision(9);

		const std::string experiment = args[1];
		if (experiment == "msd")
		{
			run_msd(walk);
		}
		else if (experiment == "tracks")
		{
			run_tracks(walk);
		}
		else
		{
			throw std::invalid_argument("unknown experiment '" + experiment + "'");
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << "motility: " << e.what() << '\n'
				  << "usage: " << args[0] << " msd|tracks [dims [persistence_time [bias]]]\n";
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
