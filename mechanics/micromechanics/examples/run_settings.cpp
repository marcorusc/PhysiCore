// Runs a sub-cellular element simulation described by a PhysiCell settings file.
//
// The file is read as described in docs/Parameters.md: domain, timestep, walls, one cell type per
// <cell_definition> with its <sub_cellular_elements>, and the positions of the cells. settings.xml and
// cells.csv next to this file are a small example with two cell types; run it from this folder, as the
// settings name the folder of cells.csv relative to where the program runs.
//
// Prints over time the number of cells and sub-agents and, per cell type, the mean speed and pressure
// of its cells. With an output folder, also writes the sub-agents as VTK files for ParaView.
//
// Usage: mechanics.micromechanics.run_settings <settings.xml> [end_time [output folder]]

#include <cmath>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include <micromechanics/environment.h>
#include <micromechanics/vtk_serializer.h>

using namespace physicore;
using namespace physicore::mechanics::micromechanics;

int main(int argc, char** argv)
{
	const std::span<char*> args(argv, static_cast<std::size_t>(argc));

	try
	{
		if (args.size() < 2 || args.size() > 4)
		{
			throw std::invalid_argument("expected a settings file and at most 2 further arguments");
		}
		const real_t end_time = args.size() > 2 ? std::stod(args[2]) : 10.0;
		if (end_time <= 0.0)
		{
			throw std::invalid_argument("the end time must be positive");
		}

		auto env = environment::create_from_config(args[1]);
		if (args.size() > 3)
		{
			std::filesystem::create_directories(args[3]);
			env->serializer = make_vtk_agents_serializer(args[3]);
		}

		std::cout << "# settings=" << args[1] << " dims=" << env->mesh.dims << " timestep=" << env->mechanics_timestep
				  << " cell_types=" << env->cell_types.size() << " cells=" << env->cells->size()
				  << " sub_agents=" << env->agents->size() << '\n';
		std::cout << "time";
		for (const cell_type& type : env->cell_types)
		{
			std::cout << ",cells_" << type.name << ",mean_speed_" << type.name << ",mean_pressure_" << type.name;
		}
		std::cout << '\n' << std::setprecision(9);

		const index_t snapshots = 50;
		const auto steps_per_snapshot = static_cast<index_t>(
			std::max(1.0, std::round(end_time / (static_cast<real_t>(snapshots) * env->mechanics_timestep))));

		for (index_t snapshot = 0; snapshot <= snapshots; ++snapshot)
		{
			if (snapshot > 0)
			{
				env->solver->solve(*env, steps_per_snapshot);
			}
			const real_t time = static_cast<real_t>(snapshot * steps_per_snapshot) * env->mechanics_timestep;
			env->serialize_state(time);

			std::vector<index_t> count(env->cell_types.size(), 0);
			std::vector<real_t> speed(env->cell_types.size(), 0.0);
			std::vector<real_t> pressure(env->cell_types.size(), 0.0);
			for (index_t c = 0; c < env->cells->size(); ++c)
			{
				auto* cell = env->cells->get_agent_at(c);
				const index_t type = cell->type_index();
				real_t squared = 0.0;
				for (const real_t v : cell->velocity())
				{
					squared += v * v;
				}
				++count[type];
				speed[type] += std::sqrt(squared);
				pressure[type] += cell->pressure();
			}

			std::cout << time;
			for (std::size_t t = 0; t < count.size(); ++t)
			{
				const real_t cells = std::max<real_t>(1.0, static_cast<real_t>(count[t]));
				std::cout << ',' << count[t] << ',' << speed[t] / cells << ',' << pressure[t] / cells;
			}
			std::cout << '\n';
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << "run_settings: " << e.what() << '\n'
				  << "usage: " << args[0] << " <settings.xml> [end_time [output folder]]\n";
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
