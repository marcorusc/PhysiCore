#include "config_reader.h"

#include <cmath>
#include <cstdint>
#include <numbers>
#include <pugixml.hpp>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include "micromechanics/cell_parameters.h"

namespace physicore::mechanics::micromechanics {

namespace {

// true for a number between 0 and 1; false also for a value that is not a number
bool is_fraction(real_t value) { return value >= 0.0 && value <= 1.0; }

pugi::xml_node required_child(const pugi::xml_node& parent, const char* name)
{
	const pugi::xml_node child = parent.child(name);
	if (!child)
	{
		throw std::runtime_error(std::string("Required XML tag <") + name + "> not found under <" + parent.name()
								 + ">");
	}
	return child;
}

real_t required_real(const pugi::xml_node& parent, const char* name)
{
	return static_cast<real_t>(required_child(parent, name).text().as_double());
}

real_t optional_real(const pugi::xml_node& parent, const char* name, real_t fallback)
{
	const pugi::xml_node child = parent.child(name);
	return child ? static_cast<real_t>(child.text().as_double()) : fallback;
}

void parse_domain(const pugi::xml_node& domain, simulation_config& config)
{
	config.dims = domain.child("use_2D").text().as_bool(false) ? 2 : 3;

	config.bounding_box_mins = { static_cast<sindex_t>(required_real(domain, "x_min")),
								 static_cast<sindex_t>(required_real(domain, "y_min")),
								 static_cast<sindex_t>(required_real(domain, "z_min")) };
	config.bounding_box_maxs = { static_cast<sindex_t>(required_real(domain, "x_max")),
								 static_cast<sindex_t>(required_real(domain, "y_max")),
								 static_cast<sindex_t>(required_real(domain, "z_max")) };

	const std::array<real_t, 3> spacing = { required_real(domain, "dx"), required_real(domain, "dy"),
											required_real(domain, "dz") };
	for (index_t d = 0; d < 3; ++d)
	{
		if (d < config.dims && config.bounding_box_maxs.at(d) <= config.bounding_box_mins.at(d))
		{
			throw std::runtime_error("<domain>: the upper bound must be above the lower bound along every axis");
		}
		if (!(spacing.at(d) >= 1.0))
		{
			throw std::runtime_error("<domain>: dx, dy and dz must be at least 1");
		}
		config.voxel_shape.at(d) = static_cast<index_t>(spacing.at(d));
	}
}

void parse_micromechanics(const pugi::xml_node& node, simulation_parameters& params)
{
	params.steepness = optional_real(node, "steepness", params.steepness);
	params.cutoff_factor = optional_real(node, "cutoff_factor", params.cutoff_factor);
	params.thermal_energy = optional_real(node, "thermal_energy", params.thermal_energy);
	if (const pugi::xml_node seed = node.child("random_seed"); seed)
	{
		params.random_seed = seed.text().as_ullong();
	}

	if (!(params.steepness > 0.0) || !(params.cutoff_factor > 0.0) || !(params.thermal_energy >= 0.0))
	{
		throw std::runtime_error("<micromechanics>: steepness and cutoff_factor must be positive and "
								 "thermal_energy must not be negative");
	}
}

// The elements of a cell type, from the properties of the whole cell or given directly
cell_type parse_elements(const pugi::xml_node& elements_node, const std::string& name, real_t radius,
						 const simulation_config& config)
{
	const auto elements = static_cast<std::int64_t>(std::llround(required_real(elements_node, "elements")));
	if (elements < 1)
	{
		throw std::runtime_error("Cell type '" + name + "': <elements> must be at least 1");
	}

	cell_type type;
	if (const pugi::xml_node element = elements_node.child("element"); element)
	{
		type.equilibrium_distance = required_real(element, "equilibrium_distance");
		type.intra_well_depth = required_real(element, "intra_well_depth");
		type.inter_well_depth = required_real(element, "inter_well_depth");
		type.drag_coefficient = required_real(element, "drag_coefficient");
		if (!(type.equilibrium_distance > 0.0) || !(type.intra_well_depth > 0.0) || !(type.inter_well_depth >= 0.0)
			|| !(type.drag_coefficient > 0.0))
		{
			throw std::runtime_error("Cell type '" + name
									 + "': the values under <element> must be positive "
									   "(inter_well_depth may be zero)");
		}
	}
	else
	{
		if (config.dims != 3)
		{
			throw std::runtime_error(
				"Cell type '" + name
				+ "': <stiffness> and <damping> of the whole cell can only be turned into elements in three "
				  "dimensions; in a 2D simulation give the elements directly under <element>");
		}

		cell_parameters cell;
		cell.radius = radius;
		cell.stiffness = required_real(elements_node, "stiffness");
		cell.damping = required_real(elements_node, "damping");
		cell.adhesion_ratio = optional_real(elements_node, "adhesion_ratio", cell.adhesion_ratio);
		try
		{
			type = scale_to_elements(cell, static_cast<index_t>(elements), config.params);
		}
		catch (const std::invalid_argument& e)
		{
			throw std::runtime_error("Cell type '" + name + "': " + e.what());
		}
	}

	type.name = name;
	type.elements_per_cell = static_cast<index_t>(elements);
	return type;
}

void parse_motility(const pugi::xml_node& motility, cell_type& type)
{
	type.migration_speed = required_real(motility, "speed");
	type.persistence_time = required_real(motility, "persistence_time");
	type.migration_bias = required_real(motility, "migration_bias");
	type.is_motile = motility.child("options").child("enabled").text().as_bool(false);

	const bool valid = type.persistence_time > 0.0 && is_fraction(type.migration_bias);
	if (type.is_motile && !valid)
	{
		throw std::runtime_error("Cell type '" + type.name
								 + "': a motile cell needs a positive persistence_time and a migration_bias "
								   "between 0 and 1");
	}
}

void parse_affinities(const pugi::xml_node& mechanics, cell_type& type,
					  const std::unordered_map<std::string, index_t>& index_of_name)
{
	// as in PhysiCell, a type that is not listed has the affinity 1
	type.adhesion_affinities.assign(index_of_name.size(), 1.0);

	const pugi::xml_node affinities = mechanics.child("cell_adhesion_affinities");
	for (pugi::xml_node affinity = affinities.child("cell_adhesion_affinity"); affinity;
		 affinity = affinity.next_sibling("cell_adhesion_affinity"))
	{
		const std::string other = affinity.attribute("name").as_string();
		const auto found = index_of_name.find(other);
		if (found == index_of_name.end())
		{
			throw std::runtime_error("Cell type '" + type.name + "': unknown cell type '" + other
									 + "' in <cell_adhesion_affinity>");
		}
		const auto value = static_cast<real_t>(affinity.text().as_double());
		if (!is_fraction(value))
		{
			throw std::runtime_error("Cell type '" + type.name + "': the affinity to '" + other
									 + "' must be between 0 and 1");
		}
		type.adhesion_affinities[found->second] = value;
	}
}

void parse_cell_definitions(const pugi::xml_node& definitions, simulation_config& config)
{
	// first the names and IDs, so that affinities can refer to types defined further down
	std::unordered_map<std::string, index_t> index_of_name;
	index_t count = 0;
	for (pugi::xml_node definition = definitions.child("cell_definition"); definition;
		 definition = definition.next_sibling("cell_definition"))
	{
		const std::string name = definition.attribute("name").as_string();
		const auto id = static_cast<index_t>(definition.attribute("ID").as_ullong(count));
		if (name.empty() || !index_of_name.emplace(name, id).second)
		{
			throw std::runtime_error("<cell_definition>: every cell type needs a name of its own");
		}
		++count;
	}
	if (count == 0)
	{
		throw std::runtime_error("<cell_definitions>: at least one <cell_definition> is needed");
	}

	// the IDs must be 0, 1, 2, ... : they are the indices of the cell types
	std::vector<bool> seen(count, false);
	for (const auto& [name, id] : index_of_name)
	{
		if (id >= count || seen[id])
		{
			throw std::runtime_error("<cell_definition>: the IDs of the cell types must be 0, 1, 2, ... without "
									 "gaps (cell type '"
									 + name + "')");
		}
		seen[id] = true;
	}

	config.cell_types.assign(count, cell_type {});

	for (pugi::xml_node definition = definitions.child("cell_definition"); definition;
		 definition = definition.next_sibling("cell_definition"))
	{
		const std::string name = definition.attribute("name").as_string();
		const index_t id = index_of_name.at(name);

		const pugi::xml_node phenotype = required_child(definition, "phenotype");
		const real_t volume = required_real(required_child(phenotype, "volume"), "total");
		if (!(volume > 0.0))
		{
			throw std::runtime_error("Cell type '" + name + "': the total volume must be positive");
		}
		// the cell as a sphere of that volume
		const real_t radius = std::cbrt(volume / (4.0 / 3.0 * std::numbers::pi));

		const pugi::xml_node mechanics = required_child(phenotype, "mechanics");
		const pugi::xml_node elements = mechanics.child("sub_cellular_elements");
		if (!elements)
		{
			throw std::runtime_error("Cell type '" + name
									 + "': <sub_cellular_elements> not found under <mechanics>; the sub-cellular "
									   "element model needs it for every cell type");
		}

		cell_type type = parse_elements(elements, name, radius, config);
		parse_affinities(mechanics, type, index_of_name);
		if (const pugi::xml_node motility = phenotype.child("motility"); motility)
		{
			parse_motility(motility, type);
		}
		config.cell_types[id] = std::move(type);
	}
}

} // namespace

simulation_config parse_simulation_config(const std::filesystem::path& config_file)
{
	pugi::xml_document doc;
	if (const pugi::xml_parse_result result = doc.load_file(config_file.string().c_str()); !result)
	{
		throw std::runtime_error("Failed to parse XML file: " + config_file.string() + " - " + result.description());
	}

	const pugi::xml_node root = doc.child("PhysiCell_settings");
	if (!root)
	{
		throw std::runtime_error("Root <PhysiCell_settings> tag not found in " + config_file.string());
	}

	simulation_config config;

	parse_domain(required_child(root, "domain"), config);

	config.timestep = required_real(required_child(root, "overall"), "dt_mechanics");
	if (!(config.timestep > 0.0))
	{
		throw std::runtime_error("<overall>: dt_mechanics must be positive");
	}

	config.params.domain_walls = root.child("options").child("virtual_wall_at_domain_edge").text().as_bool(false);
	if (const pugi::xml_node micromechanics = root.child("micromechanics"); micromechanics)
	{
		parse_micromechanics(micromechanics, config.params);
	}

	parse_cell_definitions(required_child(root, "cell_definitions"), config);

	if (const pugi::xml_node positions = root.child("initial_conditions").child("cell_positions");
		positions && positions.attribute("enabled").as_bool(true))
	{
		const std::filesystem::path folder = required_child(positions, "folder").text().as_string();
		const std::string filename = required_child(positions, "filename").text().as_string();
		config.cell_positions_file = folder / filename;
	}

	return config;
}

} // namespace physicore::mechanics::micromechanics
