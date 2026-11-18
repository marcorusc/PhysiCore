#include "micromechanics/environment_builder.h"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

#include "config_reader.h"
#include "micromechanics/solver_registry.h"
#include "micromechanics/vtk_serializer.h"

namespace physicore::mechanics::micromechanics {

real_t largest_stable_timestep(const cell_type& type, const simulation_parameters& params, index_t dims)
{
	const real_t well_depth = std::max(type.intra_well_depth, type.inter_well_depth);
	const real_t stiffness = 8.0 * well_depth * params.steepness * params.steepness
							 / (type.equilibrium_distance * type.equilibrium_distance);

	// largest eigenvalue of the stiffness matrix of a close packing, in units of the pair stiffness
	const real_t largest_eigenvalue = 2.0 * (static_cast<real_t>(dims) + 1.0); // 4, 6, 8

	return 2.0 * type.drag_coefficient / (largest_eigenvalue * stiffness);
}

void environment_builder::set_time_step(real_t time_step) { timestep = time_step; }

void environment_builder::resize(index_t dims, std::array<sindex_t, 3> bounding_box_mins,
								 std::array<sindex_t, 3> bounding_box_maxs, std::array<index_t, 3> voxel_shape)
{
	mesh = cartesian_mesh(dims, bounding_box_mins, bounding_box_maxs, voxel_shape);
}

void environment_builder::set_parameters(const simulation_parameters& parameters) { params = parameters; }

void environment_builder::set_domain_walls(bool value) { params.domain_walls = value; }

index_t environment_builder::add_cell_type(cell_type type)
{
	cell_types.push_back(std::move(type));
	return cell_types.size() - 1;
}

void environment_builder::set_solver(std::string name) { solver_name = std::move(name); }

void environment_builder::set_vtk_output(std::string output_dir) { vtk_output = std::move(output_dir); }

std::unique_ptr<environment> environment_builder::build()
{
	if (!mesh)
	{
		throw std::invalid_argument("The domain is not set: call resize() before build()");
	}
	if (!(timestep > 0.0))
	{
		throw std::invalid_argument("The timestep must be positive: call set_time_step() before build()");
	}

	auto env = std::make_unique<environment>(*mesh, timestep);
	env->params = params;
	if (!cell_types.empty())
	{
		env->cell_types = cell_types;
	}

	// checks the values of every type
	const interaction_table table(env->cell_types, env->params);

	for (const cell_type& type : env->cell_types)
	{
		const real_t limit = largest_stable_timestep(type, env->params, mesh->dims);
		if (timestep >= limit)
		{
			std::ostringstream message;
			message << "The timestep " << timestep << " is too large for cell type '" << type.name
					<< "': it must stay below " << limit << ", and well below it to be accurate";
			throw std::invalid_argument(message.str());
		}
	}

	env->solver = solver_registry::instance().get(solver_name);
	if (!env->solver)
	{
		throw std::invalid_argument("The solver '" + solver_name + "' is not registered");
	}
	env->solver->initialize(*env);

	if (!vtk_output.empty())
	{
		env->serializer = make_vtk_agents_serializer(vtk_output);
	}

	return env;
}

// ------------------------------------------------------------------------------- from a settings file

namespace {

// The whole of text as a number. (std::from_chars would do, but not every standard library has it
// for floating-point numbers.)
bool parse_real(const std::string& text, real_t& value)
{
	if (text.empty())
	{
		return false;
	}
	char* end = nullptr;
	value = std::strtod(text.c_str(), &end);
	return end == text.c_str() + text.size();
}

// The whole of text as a count: digits only
bool parse_index(const std::string& text, index_t& value)
{
	if (text.empty() || text.find_first_not_of("0123456789") != std::string::npos)
	{
		return false;
	}
	value = std::strtoull(text.c_str(), nullptr, 10);
	return true;
}

// One cell per line: x, y, z, type. The type is the ID of a cell type or its name. A first line that
// does not start with a number is taken as a header. Further columns are ignored.
void load_cell_positions(environment& env, const std::filesystem::path& file)
{
	std::ifstream input(file);
	if (!input)
	{
		throw std::runtime_error("Cannot open the cell positions file " + file.string());
	}

	const auto fail = [&file](index_t line_number, const std::string& what) {
		throw std::runtime_error(file.string() + ", line " + std::to_string(line_number) + ": " + what);
	};

	std::string line;
	index_t line_number = 0;
	while (std::getline(input, line))
	{
		++line_number;
		if (!line.empty() && line.back() == '\r')
		{
			line.pop_back();
		}
		if (line.find_first_not_of(" \t") == std::string::npos)
		{
			continue;
		}

		std::vector<std::string> fields;
		std::stringstream stream(line);
		for (std::string field; std::getline(stream, field, ',');)
		{
			const auto first = field.find_first_not_of(" \t");
			const auto last = field.find_last_not_of(" \t");
			fields.push_back(first == std::string::npos ? std::string {} : field.substr(first, last - first + 1));
		}

		std::array<real_t, 3> position {};
		bool numeric = fields.size() >= 4;
		for (std::size_t d = 0; d < 3 && numeric; ++d)
		{
			numeric = parse_real(fields[d], position.at(d));
		}
		if (!numeric)
		{
			if (line_number == 1)
			{
				continue; // the header
			}
			fail(line_number, "expected x, y, z, type");
		}

		// the type: an ID, or else a name
		index_t type_index = env.cell_types.size();
		const std::string& type_field = fields[3];
		if (!parse_index(type_field, type_index))
		{
			const auto found = std::ranges::find(env.cell_types, type_field, &cell_type::name);
			type_index = static_cast<index_t>(found - env.cell_types.begin());
		}
		if (type_index >= env.cell_types.size())
		{
			fail(line_number, "unknown cell type '" + type_field + "'");
		}

		env.create_cell_of_type(type_index, std::span<const real_t>(position.data(), env.mesh.dims));
	}
}

} // namespace

std::unique_ptr<environment> environment::create_from_config(const std::filesystem::path& config_file)
{
	const simulation_config config = parse_simulation_config(config_file);

	environment_builder builder;
	builder.set_time_step(config.timestep);
	builder.resize(config.dims, config.bounding_box_mins, config.bounding_box_maxs, config.voxel_shape);
	builder.set_parameters(config.params);
	for (const cell_type& type : config.cell_types)
	{
		builder.add_cell_type(type);
	}

	std::unique_ptr<environment> env;
	try
	{
		env = builder.build();
		if (config.cell_positions_file)
		{
			load_cell_positions(*env, *config.cell_positions_file);
			// brings the cells in line with their sub-agents and checks that all are inside the domain
			env->solver->initialize(*env);
		}
	}
	catch (const std::invalid_argument& e)
	{
		throw std::runtime_error(config_file.string() + ": " + e.what());
	}

	return env;
}

} // namespace physicore::mechanics::micromechanics
