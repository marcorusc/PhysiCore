// Two types of cells that prefer their own kind, in two or three dimensions.
//
// Cells of 19 elements each start as a block (a 6 x 6 sheet in 2D, a 4 x 4 x 4 block in 3D) in which
// the two types are mixed at random. Cells of the same type stick to each other with the affinity 1,
// cells of different types with a lower affinity. Something must let the cells change neighbours:
// crawling in random directions (migration speed), thermal noise (thermal energy), or both. With a low
// affinity between the types the cells then sort: contacts between unlike cells become rare
// (differential adhesion, Steinberg 1963).
//
// Prints over time two measures of mixing and how many elements have strayed from their cell (a sign
// that the noise is tearing cells apart):
//
// - nearest neighbours: for every cell, the fraction of its nearest cells (6 in 2D, 12 in 3D, by
//   the distance between centres) that are of the other type, averaged over all cells. It does not
//   ask whether cells touch, only who is next to whom, and is the measure of sorting used here.
// - contacts: pairs of cells with two elements closer than 1.3 equilibrium distances. Unlike cells
//   that do not stick stop touching long before they change places, so this falls much faster
//   and says little about where the cells are. With an output folder,
//   also writes the sub-agents as VTK files for ParaView.
//
// Usage: mechanics.micromechanics.cell_sorting [affinity_between_types [end_time [output folder]]]
//                                              [dims=2] [kT=0] [speed=0.3] [stiffness=1] [box=...]
//
//   dims       2 or 3
//   kT         thermal energy; compare with the well depths, 1 inside a cell and 0.5 between cells
//   speed      migration speed of all cells; 0 switches the crawling off
//   box        half the side of the walled box, a whole number; the default (22 in 2D, 11 in 3D) leaves
//              the cells much room. A smaller box packs them: the block of cells is then created
//              with a smaller spacing so that it fits.
//   stiffness  intra-cellular well depth of both types. The inter-cellular well stays at 0.5, so this
//              changes how easily a cell deforms without changing how strongly cells stick.

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
#include <vector>

#include <micromechanics/counter_random.h>
#include <micromechanics/environment.h>
#include <micromechanics/solver_registry.h>
#include <micromechanics/vtk_serializer.h>

using namespace physicore;
using namespace physicore::mechanics::micromechanics;

namespace {

constexpr index_t elements_per_cell = 19;
// timestep for the stiffness 1; it is divided by the stiffness, because a stiffer well is a faster
// spring and forward Euler needs to resolve it. A snapshot is taken every two time units.
constexpr real_t base_timestep = 0.005;
constexpr real_t snapshot_interval = 2.0;
// two elements of different cells closer than this count as a contact between the cells
constexpr real_t contact_distance = 1.3;
// an element further than this from the centre of its cell has strayed (a cell reaches about 2)
constexpr real_t stray_distance = 4.0;

struct options
{
	real_t affinity = 0.1;
	real_t end_time = 200.0;
	std::string output;
	index_t dims = 2;
	real_t thermal_energy = 0.0;
	real_t speed = 0.3;
	real_t stiffness = 1.0;
	sindex_t half_box = 0; // 0: the default for the dimension
};

options parse(std::span<char*> args)
{
	options opts;
	index_t positional = 0;
	for (std::size_t k = 1; k < args.size(); ++k)
	{
		const std::string arg = args[k];
		const auto equals = arg.find('=');
		if (equals == std::string::npos)
		{
			switch (positional++)
			{
				case 0:
					opts.affinity = std::stod(arg);
					break;
				case 1:
					opts.end_time = std::stod(arg);
					break;
				case 2:
					opts.output = arg;
					break;
				default:
					throw std::invalid_argument("too many arguments");
			}
			continue;
		}

		const std::string key = arg.substr(0, equals);
		const std::string value = arg.substr(equals + 1);
		if (key == "dims")
			opts.dims = std::stoull(value);
		else if (key == "kT")
			opts.thermal_energy = std::stod(value);
		else if (key == "speed")
			opts.speed = std::stod(value);
		else if (key == "stiffness")
			opts.stiffness = std::stod(value);
		else if (key == "box")
			opts.half_box = std::stoll(value);
		else
			throw std::invalid_argument("unknown option '" + key + "'");
	}

	if (opts.affinity < 0.0 || opts.affinity > 1.0 || opts.end_time <= 0.0 || (opts.dims != 2 && opts.dims != 3)
		|| opts.thermal_energy < 0.0 || opts.speed < 0.0 || opts.stiffness <= 0.0)
	{
		throw std::invalid_argument("affinity in [0, 1], positive end time and stiffness, dims 2 or 3, "
									"kT and speed not negative");
	}
	return opts;
}

struct statistics
{
	index_t contacts = 0;
	index_t unlike = 0;
	real_t unlike_nearest = 0.0;
	index_t strays = 0;
};

// Who is next to whom: for every cell, the fraction of its nearest cells (by the distance between
// centres) that are of the other type, averaged over all cells
real_t unlike_among_nearest(environment& env)
{
	const index_t dims = env.mesh.dims;
	real_t average = 0.0;

	const index_t cells_count = env.cells->size();
	const index_t nearest = std::min<index_t>(dims == 2 ? 6 : 12, cells_count - 1);
	std::vector<std::pair<real_t, index_t>> others;
	for (index_t a = 0; a < cells_count; ++a)
	{
		auto* cell = env.cells->get_agent_at(a);
		others.clear();
		for (index_t b = 0; b < cells_count; ++b)
		{
			if (a == b)
			{
				continue;
			}
			auto* other = env.cells->get_agent_at(b);
			real_t squared = 0.0;
			for (index_t d = 0; d < dims; ++d)
			{
				squared += std::pow(cell->position()[d] - other->position()[d], 2);
			}
			others.emplace_back(squared, other->type_index());
		}
		std::ranges::partial_sort(others, others.begin() + static_cast<std::ptrdiff_t>(nearest));
		index_t unlike = 0;
		for (index_t k = 0; k < nearest; ++k)
		{
			unlike += others[k].second != cell->type_index() ? 1 : 0;
		}
		average += static_cast<real_t>(unlike) / static_cast<real_t>(nearest) / static_cast<real_t>(cells_count);
	}
	return average;
}

statistics measure(environment& env)
{
	const index_t dims = env.mesh.dims;
	const index_t count = env.agents->size();

	// copied out once: the loop over pairs below is the expensive part
	std::vector<real_t> positions(count * dims);
	std::vector<index_t> cell_of(count);
	for (index_t i = 0; i < count; ++i)
	{
		auto* element = env.agents->get_agent_at(i);
		std::ranges::copy(element->position(), positions.begin() + static_cast<std::ptrdiff_t>(i * dims));
		cell_of[i] = element->cell_index();
	}

	statistics result;
	std::set<std::pair<index_t, index_t>> touching;
	for (index_t i = 0; i < count; ++i)
	{
		const auto centre = env.cells->get_agent_at(cell_of[i])->position();
		real_t from_centre = 0.0;
		for (index_t d = 0; d < dims; ++d)
		{
			from_centre += std::pow(positions[i * dims + d] - centre[d], 2);
		}
		if (from_centre > stray_distance * stray_distance)
		{
			++result.strays;
		}

		for (index_t j = i + 1; j < count; ++j)
		{
			if (cell_of[i] == cell_of[j])
			{
				continue;
			}
			real_t squared = 0.0;
			for (index_t d = 0; d < dims; ++d)
			{
				squared += std::pow(positions[i * dims + d] - positions[j * dims + d], 2);
			}
			if (squared < contact_distance * contact_distance)
			{
				touching.insert(std::minmax(cell_of[i], cell_of[j]));
			}
		}
	}

	result.unlike_nearest = unlike_among_nearest(env);

	result.contacts = touching.size();
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
		const options opts = parse(args);
		const real_t timestep = base_timestep / std::max(1.0, opts.stiffness);
		const auto steps_per_snapshot = static_cast<index_t>(std::llround(snapshot_interval / timestep));

		// the block of cells: 6 x 6 cells 5 apart in 2D (a disc of 19 elements reaches 2 from its
		// centre), 4 x 4 x 4 cells 3.6 apart in 3D (a ball of 19 elements reaches 1.4)
		const index_t cells_per_side = opts.dims == 2 ? 6 : 4;
		const sindex_t default_half_box = opts.dims == 2 ? 22 : 11;
		const sindex_t half_box = opts.half_box > 0 ? opts.half_box : default_half_box;
		// in a small box the block is created tighter, so that its outermost elements start about half
		// an equilibrium distance inside the walls
		const real_t reach = opts.dims == 2 ? 2.0 : 1.5;
		const real_t cell_spacing =
			std::min(opts.dims == 2 ? 5.0 : 3.6,
					 2.0 * (static_cast<real_t>(half_box) - reach - 0.5) / static_cast<real_t>(cells_per_side - 1));

		// a box with walls, so that the cells stay together
		environment env(
			{ opts.dims, { -half_box, -half_box, -half_box }, { half_box, half_box, half_box }, { 1, 1, 1 } },
			timestep);
		env.params.domain_walls = true;
		env.params.random_seed = 2026;
		env.params.thermal_energy = opts.thermal_energy;

		cell_type first;
		first.name = "first";
		first.elements_per_cell = elements_per_cell;
		first.intra_well_depth = opts.stiffness;
		first.is_motile = opts.speed > 0.0;
		first.migration_speed = opts.speed;
		first.persistence_time = 2.0;
		first.adhesion_affinities = { 1.0, opts.affinity };

		cell_type second = first;
		second.name = "second";
		second.adhesion_affinities = { opts.affinity, 1.0 };

		env.cell_types[0] = first;
		env.add_cell_type(second);

		// each cell is of a type drawn at random with equal chances
		const index_t layers = opts.dims == 2 ? 1 : cells_per_side;
		const real_t offset = 0.5 * static_cast<real_t>(cells_per_side - 1);
		index_t created = 0;
		for (index_t layer = 0; layer < layers; ++layer)
		{
			for (index_t row = 0; row < cells_per_side; ++row)
			{
				for (index_t column = 0; column < cells_per_side; ++column)
				{
					const std::array<real_t, 3> centre = { cell_spacing * (static_cast<real_t>(column) - offset),
														   cell_spacing * (static_cast<real_t>(row) - offset),
														   cell_spacing * (static_cast<real_t>(layer) - offset) };
					const index_t type = random_bits(env.params.random_seed, 0, created++, 0, 7) % 2;
					env.create_cell_of_type(type, std::span<const real_t>(centre.data(), opts.dims));
				}
			}
		}

		env.solver = solver_registry::instance().get("openmp_solver");
		if (!env.solver)
		{
			throw std::runtime_error("the openmp_solver is not registered");
		}
		env.solver->initialize(env);

		if (!opts.output.empty())
		{
			std::filesystem::create_directories(opts.output);
			env.serializer = make_vtk_agents_serializer(opts.output);
		}

		index_t of_second_type = 0;
		for (index_t c = 0; c < env.cells->size(); ++c)
		{
			of_second_type += env.cells->get_agent_at(c)->type_index();
		}

		std::cout << std::setprecision(9);
		std::cout << "# affinity_between_types=" << opts.affinity << " cells=" << env.cells->size()
				  << " of_second_type=" << of_second_type << " elements_per_cell=" << elements_per_cell
				  << " dims=" << opts.dims << " thermal_energy=" << opts.thermal_energy
				  << " migration_speed=" << opts.speed << " persistence_time=" << first.persistence_time
				  << " intra_well_depth=" << opts.stiffness << " inter_well_depth=" << first.inter_well_depth
				  << " half_box=" << half_box << " cell_spacing=" << cell_spacing << " timestep=" << timestep << '\n';
		std::cout << "time,unlike_nearest_fraction,contacts,unlike_contacts,"
					 "unlike_fraction,stray_elements\n";

		const auto snapshots =
			static_cast<index_t>(std::llround(opts.end_time / (timestep * static_cast<real_t>(steps_per_snapshot))));
		for (index_t snapshot = 0; snapshot <= snapshots; ++snapshot)
		{
			if (snapshot > 0)
			{
				env.solver->solve(env, steps_per_snapshot);
			}
			const real_t time = static_cast<real_t>(snapshot * steps_per_snapshot) * timestep;
			env.serialize_state(time);

			const statistics now = measure(env);
			const auto fraction = [](index_t part, index_t whole) {
				return whole == 0 ? 0.0 : static_cast<real_t>(part) / static_cast<real_t>(whole);
			};
			std::cout << time << ',' << now.unlike_nearest << ',' << now.contacts << ',' << now.unlike << ','
					  << fraction(now.unlike, now.contacts) << ',' << now.strays << std::endl;
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << "cell_sorting: " << e.what() << '\n'
				  << "usage: " << args[0] << " [affinity_between_types [end_time [output folder]]]"
				  << " [dims=2] [kT=0] [speed=0.3] [stiffness=1]\n";
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
