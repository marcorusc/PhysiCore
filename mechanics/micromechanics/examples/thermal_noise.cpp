// Two small experiments on the thermal noise, printed as CSV.
//
//   diffusion [dims [diffusion_coefficient]]
//       Free elements that do not interact. Prints the mean squared displacement over time, which
//       should grow as 2 * dims * D * time.
//
//   pair [thermal_energy_over_stiffness]
//       Bound pairs of elements in thermal equilibrium. Prints sampled pair distances, which should
//       follow the Boltzmann distribution of the pair potential at kT = D * drag_coefficient.
//
// Usage: mechanics.micromechanics.thermal_noise diffusion|pair [...]

#include <cmath>
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <span>
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

void run_diffusion(index_t dims, real_t diffusion)
{
	const index_t elements = 4000;
	const real_t timestep = 0.01;
	const index_t snapshots = 40;
	const index_t steps_per_snapshot = 10;
	const real_t separation = 10.0;

	environment env = make_environment(dims, timestep);
	env.params.thermal_energy = diffusion;
	// elements without a cell use the inter-cellular potential; switching it off makes them free
	env.cell_types[0].inter_well_depth = 0.0;
	for (index_t i = 0; i < elements; ++i)
	{
		env.agents->create()->position()[0] = separation * static_cast<real_t>(i);
	}
	attach_solver(env);

	std::cout << "# experiment=diffusion dims=" << dims << " diffusion_coefficient=" << diffusion
			  << " elements=" << elements << " timestep=" << timestep << '\n';
	std::cout << "time,mean_squared_displacement\n";

	for (index_t snapshot = 0; snapshot <= snapshots; ++snapshot)
	{
		if (snapshot > 0)
		{
			env.solver->solve(env, steps_per_snapshot);
		}

		real_t sum = 0.0;
		for (index_t i = 0; i < elements; ++i)
		{
			const auto position = env.agents->get_agent_at(i)->position();
			for (index_t d = 0; d < dims; ++d)
			{
				const real_t start = d == 0 ? separation * static_cast<real_t>(i) : 0.0;
				sum += std::pow(position[d] - start, 2);
			}
		}
		std::cout << static_cast<real_t>(snapshot * steps_per_snapshot) * timestep << ','
				  << sum / static_cast<real_t>(elements) << '\n';
	}
}

void run_pair(real_t thermal_energy_over_stiffness)
{
	const index_t pairs = 400;
	const index_t samples = 20;

	const simulation_parameters params;
	const cell_type type;
	const real_t stiffness = 8.0 * type.intra_well_depth * params.steepness * params.steepness
							 / (type.equilibrium_distance * type.equilibrium_distance);
	const real_t thermal_energy = thermal_energy_over_stiffness * stiffness;

	// the pair distance relaxes in drag / (2 k); the timestep resolves that 25 times
	const real_t timestep = 0.02 * type.drag_coefficient / stiffness;

	environment env = make_environment(3, timestep);
	env.params.thermal_energy = thermal_energy;
	env.params.random_seed = 2026;
	env.create_cell();
	for (index_t p = 0; p < pairs; ++p)
	{
		for (const real_t y : { 0.0, type.equilibrium_distance })
		{
			auto position = env.create_agent(0)->position();
			position[0] = 10.0 * static_cast<real_t>(p);
			position[1] = y;
		}
	}
	attach_solver(env);

	std::cout << "# experiment=pair thermal_energy=" << thermal_energy << " stiffness=" << stiffness
			  << " well_depth=" << type.intra_well_depth << " equilibrium_distance=" << type.equilibrium_distance
			  << " steepness=" << params.steepness << " cutoff_factor=" << params.cutoff_factor
			  << " drag_coefficient=" << type.drag_coefficient << " timestep=" << timestep << '\n';
	std::cout << "distance\n";
	// nine digits are plenty for a histogram and keep the output small
	std::cout << std::setprecision(9);

	env.solver->solve(env, 500); // twenty relaxation times to forget the start
	for (index_t sample = 0; sample < samples; ++sample)
	{
		env.solver->solve(env, 50); // two relaxation times between samples
		for (index_t p = 0; p < pairs; ++p)
		{
			const auto a = env.agents->get_agent_at(2 * p)->position();
			const auto b = env.agents->get_agent_at(2 * p + 1)->position();
			std::cout << std::sqrt(std::pow(a[0] - b[0], 2) + std::pow(a[1] - b[1], 2) + std::pow(a[2] - b[2], 2))
					  << '\n';
		}
	}
}

} // namespace

int main(int argc, char** argv)
{
	const std::span<char*> args(argv, static_cast<std::size_t>(argc));

	try
	{
		if (args.size() < 2 || args.size() > 4)
		{
			throw std::invalid_argument("expected an experiment and at most 2 further arguments");
		}

		std::cout << std::setprecision(std::numeric_limits<real_t>::max_digits10);

		const std::string experiment = args[1];
		if (experiment == "diffusion")
		{
			const index_t dims = args.size() > 2 ? std::stoull(args[2]) : 3;
			const real_t diffusion = args.size() > 3 ? std::stod(args[3]) : 0.16;
			if (dims < 1 || dims > 3 || diffusion <= 0.0)
			{
				throw std::invalid_argument("dims must be 1, 2 or 3 and the diffusion coefficient positive");
			}
			run_diffusion(dims, diffusion);
		}
		else if (experiment == "pair")
		{
			const real_t ratio = args.size() > 2 ? std::stod(args[2]) : 1e-3;
			if (ratio <= 0.0 || args.size() > 3)
			{
				throw std::invalid_argument("pair takes one positive value");
			}
			run_pair(ratio);
		}
		else
		{
			throw std::invalid_argument("unknown experiment '" + experiment + "'");
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << "thermal_noise: " << e.what() << '\n'
				  << "usage: " << args[0] << " diffusion [dims [diffusion_coefficient]]\n"
				  << "       " << args[0] << " pair [thermal_energy_over_stiffness]\n";
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
