// Two types of motile cells that prefer their own kind, in two dimensions.
//
// 36 cells of 19 elements each start as a sheet in which the two types are mixed at random. Cells of
// the same type stick to each other with the affinity 1, cells of different types with the affinity
// given on the command line. All cells crawl in random directions, which lets them change neighbours.
// With a low affinity between the types the cells sort: contacts between unlike cells become rare
// (differential adhesion, Steinberg 1963).
//
// Prints over time the number of pairs of cells in contact and the fraction of them that joins cells of
// different types. With an output folder, also writes the sub-agents as VTK files for ParaView.
//
// Usage: mechanics.micromechanics.cell_sorting [affinity_between_types [end_time [output folder]]]

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <set>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>

#include <micromechanics/counter_random.h>
#include <micromechanics/environment.h>
#include <micromechanics/solver_registry.h>
#include <micromechanics/vtk_serializer.h>

using namespace physicore;
using namespace physicore::mechanics::micromechanics;

namespace {

constexpr index_t cells_per_side = 6;
constexpr index_t elements_per_cell = 19; // a hexagonal disc: the centre and two rings
constexpr real_t cell_spacing = 5.0;	  // the disc reaches 2 from its centre
constexpr real_t timestep = 0.005;
constexpr index_t steps_per_snapshot = 400;
// two elements of different cells closer than this count as a contact between the cells
constexpr real_t contact_distance = 1.3;

struct contacts
{
	index_t all = 0;
	index_t unlike = 0;
};

contacts count_contacts(environment& env)
{
	std::set<std::pair<index_t, index_t>> touching;
	const index_t count = env.agents->size();
	for (index_t i = 0; i < count; ++i)
	{
		auto* a = env.agents->get_agent_at(i);
		for (index_t j = i + 1; j < count; ++j)
		{
			auto* b = env.agents->get_agent_at(j);
			if (a->cell_index() == b->cell_index())
			{
				continue;
			}
			const real_t distance =
				std::hypot(a->position()[0] - b->position()[0], a->position()[1] - b->position()[1]);
			if (distance < contact_distance)
			{
				touching.insert(std::minmax(a->cell_index(), b->cell_index()));
			}
		}
	}

	contacts result;
	result.all = touching.size();
	for (const auto& [first, second] : touching)
	{
		if (env.cells->get_agent_at(first)->type_index() != env.cells->get_agent_at(second)->type_index())
		{
			++result.unlike;
		}
	}
	return result;
}

} // namespace

int main(int argc, char** argv)
{
	const std::span<char*> args(argv, static_cast<std::size_t>(argc));

	try
	{
		if (args.size() > 4)
		{
			throw std::invalid_argument("expected at most 3 arguments");
		}
		const real_t affinity = args.size() > 1 ? std::stod(args[1]) : 0.1;
		const real_t end_time = args.size() > 2 ? std::stod(args[2]) : 200.0;
		if (affinity < 0.0 || affinity > 1.0 || end_time <= 0.0)
		{
			throw std::invalid_argument("the affinity must be between 0 and 1 and the end time positive");
		}

		// a box with walls, so that the cells stay together
		environment env({ 2, { -22, -22, -22 }, { 22, 22, 22 }, { 1, 1, 1 } }, timestep);
		env.params.domain_walls = true;
		env.params.random_seed = 2026;

		cell_type first;
		first.name = "first";
		first.elements_per_cell = elements_per_cell;
		first.is_motile = true;
		first.migration_speed = 0.3;
		first.persistence_time = 2.0;
		first.adhesion_affinities = { 1.0, affinity };

		cell_type second = first;
		second.name = "second";
		second.adhesion_affinities = { affinity, 1.0 };

		env.cell_types[0] = first;
		env.add_cell_type(second);

		// a sheet of cells, each of a type drawn at random with equal chances
		for (index_t row = 0; row < cells_per_side; ++row)
		{
			for (index_t column = 0; column < cells_per_side; ++column)
			{
				const std::array<real_t, 2> centre = {
					cell_spacing * (static_cast<real_t>(column) - 0.5 * (cells_per_side - 1)),
					cell_spacing * (static_cast<real_t>(row) - 0.5 * (cells_per_side - 1))
				};
				const index_t type = random_bits(env.params.random_seed, 0, row * cells_per_side + column, 0, 7) % 2;
				env.create_cell_of_type(type, centre);
			}
		}

		env.solver = solver_registry::instance().get("openmp_solver");
		if (!env.solver)
		{
			throw std::runtime_error("the openmp_solver is not registered");
		}
		env.solver->initialize(env);

		if (args.size() > 3)
		{
			const std::filesystem::path output = args[3];
			std::filesystem::create_directories(output);
			env.serializer = make_vtk_agents_serializer(output.string());
		}

		index_t of_second_type = 0;
		for (index_t c = 0; c < env.cells->size(); ++c)
		{
			of_second_type += env.cells->get_agent_at(c)->type_index();
		}

		std::cout << std::setprecision(9);
		std::cout << "# affinity_between_types=" << affinity << " cells=" << env.cells->size()
				  << " of_second_type=" << of_second_type << " elements_per_cell=" << elements_per_cell
				  << " migration_speed=" << first.migration_speed << " persistence_time=" << first.persistence_time
				  << " timestep=" << timestep << '\n';
		std::cout << "time,contacts,unlike_contacts,unlike_fraction\n";

		const auto snapshots = static_cast<index_t>(std::llround(end_time / (timestep * steps_per_snapshot)));
		for (index_t snapshot = 0; snapshot <= snapshots; ++snapshot)
		{
			if (snapshot > 0)
			{
				env.solver->solve(env, steps_per_snapshot);
			}
			const real_t time = static_cast<real_t>(snapshot * steps_per_snapshot) * timestep;
			env.serialize_state(time);

			const contacts now = count_contacts(env);
			std::cout << time << ',' << now.all << ',' << now.unlike << ','
					  << (now.all == 0 ? 0.0 : static_cast<real_t>(now.unlike) / static_cast<real_t>(now.all)) << '\n';
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << "cell_sorting: " << e.what() << '\n'
				  << "usage: " << args[0] << " [affinity_between_types [end_time [output folder]]]\n";
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
