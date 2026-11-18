#include <cmath>
#include <filesystem>
#include <fstream>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "config_reader.h"
#include "micromechanics/cell_parameters.h"
#include "micromechanics/environment.h"
#include "micromechanics/environment_builder.h"

// A simulation is described by the PhysiCell settings file, with one block per cell type for the
// sub-cellular elements. These tests check what the reader takes from the file, what it refuses, and
// that an environment made from a file is the one made by hand from the same values.

namespace physicore::mechanics::micromechanics::tests {

namespace {

// The settings of a simulation, as text that single tests change before writing it to a file
struct settings
{
	std::string domain = R"(
    <domain>
        <x_min>-40</x_min> <x_max>40</x_max>
        <y_min>-30</y_min> <y_max>30</y_max>
        <z_min>-20</z_min> <z_max>20</z_max>
        <dx>20</dx> <dy>20</dy> <dz>20</dz>
        <use_2D>false</use_2D>
    </domain>)";

	std::string overall = R"(
    <overall>
        <max_time units="min">30</max_time>
        <dt_mechanics units="min">0.0001</dt_mechanics>
    </overall>)";

	std::string options = R"(
    <options>
        <virtual_wall_at_domain_edge>true</virtual_wall_at_domain_edge>
    </options>)";

	std::string micromechanics = R"(
    <micromechanics>
        <steepness>3</steepness>
        <cutoff_factor>2.2</cutoff_factor>
        <thermal_energy>0.04</thermal_energy>
        <random_seed>77</random_seed>
    </micromechanics>)";

	// the cell of Sandersius and Newman (2008): a sphere of radius 10 has the volume 4188.79
	std::string epithelial_elements = R"(
                    <sub_cellular_elements>
                        <elements>128</elements>
                        <stiffness>5000</stiffness>
                        <damping>5000</damping>
                        <adhesion_ratio>0.4</adhesion_ratio>
                    </sub_cellular_elements>)";

	std::string epithelial_affinities = R"(
                    <cell_adhesion_affinities>
                        <cell_adhesion_affinity name="fibroblast">0.3</cell_adhesion_affinity>
                    </cell_adhesion_affinities>)";

	std::string fibroblast = R"(
        <cell_definition name="fibroblast" ID="1">
            <phenotype>
                <volume><total>523.598776</total></volume>
                <mechanics>
                    <cell_cell_adhesion_strength>0.4</cell_cell_adhesion_strength>
                    <cell_adhesion_affinities>
                        <cell_adhesion_affinity name="epithelial">0.2</cell_adhesion_affinity>
                        <cell_adhesion_affinity name="fibroblast">0.9</cell_adhesion_affinity>
                    </cell_adhesion_affinities>
                    <sub_cellular_elements>
                        <elements>19</elements>
                        <element>
                            <equilibrium_distance>1.5</equilibrium_distance>
                            <intra_well_depth>2</intra_well_depth>
                            <inter_well_depth>0.6</inter_well_depth>
                            <drag_coefficient>3</drag_coefficient>
                        </element>
                    </sub_cellular_elements>
                </mechanics>
                <motility>
                    <speed>0.7</speed>
                    <persistence_time>5</persistence_time>
                    <migration_bias>0.25</migration_bias>
                    <options><enabled>true</enabled></options>
                </motility>
            </phenotype>
        </cell_definition>)";

	std::string initial_conditions;

	std::string text() const
	{
		return std::string("<?xml version=\"1.0\"?>\n<PhysiCell_settings version=\"devel-version\">") + domain + overall
			   + options + micromechanics + R"(
    <cell_definitions>
        <cell_definition name="epithelial" ID="0">
            <phenotype>
                <cycle code="5" name="live"/>
                <volume><total>4188.790205</total><fluid_fraction>0.75</fluid_fraction></volume>
                <mechanics>)"
			   + epithelial_affinities + epithelial_elements + R"(
                </mechanics>
            </phenotype>
        </cell_definition>)"
			   + fibroblast + "\n    </cell_definitions>" + initial_conditions + "\n</PhysiCell_settings>\n";
	}
};

class ConfigReaderTest : public ::testing::Test
{
protected:
	std::filesystem::path folder;

	void SetUp() override
	{
		folder =
			std::filesystem::temp_directory_path()
			/ ("physicore_sem_config_" + std::string(::testing::UnitTest::GetInstance()->current_test_info()->name()));
		std::filesystem::remove_all(folder);
		std::filesystem::create_directories(folder);
	}

	void TearDown() override { std::filesystem::remove_all(folder); }

	std::filesystem::path write(const std::string& name, const std::string& content) const
	{
		const std::filesystem::path file = folder / name;
		std::ofstream(file) << content;
		return file;
	}

	std::filesystem::path write(const settings& s) const { return write("settings.xml", s.text()); }

	// the message of the std::runtime_error that reading the file throws, empty if it does not throw
	std::string error_of(const settings& s) const
	{
		try
		{
			parse_simulation_config(write(s));
		}
		catch (const std::runtime_error& e)
		{
			return e.what();
		}
		return {};
	}
};

void replace(std::string& text, const std::string& old_part, const std::string& new_part)
{
	const auto at = text.find(old_part);
	ASSERT_NE(at, std::string::npos) << old_part;
	text.replace(at, old_part.size(), new_part);
}

} // namespace

// --------------------------------------------------------------------------------- what is read

TEST_F(ConfigReaderTest, ReadsDomainTimestepWallsAndSharedParameters)
{
	const simulation_config config = parse_simulation_config(write(settings {}));

	EXPECT_EQ(config.dims, index_t { 3 });
	EXPECT_EQ(config.bounding_box_mins, (std::array<sindex_t, 3> { -40, -30, -20 }));
	EXPECT_EQ(config.bounding_box_maxs, (std::array<sindex_t, 3> { 40, 30, 20 }));
	EXPECT_EQ(config.voxel_shape, (std::array<index_t, 3> { 20, 20, 20 }));
	EXPECT_DOUBLE_EQ(config.timestep, 0.0001);

	EXPECT_TRUE(config.params.domain_walls);
	EXPECT_DOUBLE_EQ(config.params.steepness, 3.0);
	EXPECT_DOUBLE_EQ(config.params.cutoff_factor, 2.2);
	EXPECT_DOUBLE_EQ(config.params.thermal_energy, 0.04);
	EXPECT_EQ(config.params.random_seed, 77U);

	EXPECT_FALSE(config.cell_positions_file.has_value());
}

TEST_F(ConfigReaderTest, CellGivenByStiffnessAndDampingIsScaledToItsElements)
{
	const simulation_config config = parse_simulation_config(write(settings {}));
	ASSERT_EQ(config.cell_types.size(), 2U);
	const cell_type& epithelial = config.cell_types[0];

	// the radius is that of the sphere with the total volume; the rest is scale_to_elements with the
	// steepness of the file
	cell_parameters cell { .radius = 10.0, .stiffness = 5000.0, .damping = 5000.0 };
	cell.adhesion_ratio = 0.4;
	const cell_type expected = scale_to_elements(cell, 128, config.params);

	EXPECT_EQ(epithelial.name, "epithelial");
	EXPECT_EQ(epithelial.elements_per_cell, index_t { 128 });
	EXPECT_NEAR(epithelial.equilibrium_distance, expected.equilibrium_distance, 1e-8);
	EXPECT_NEAR(epithelial.intra_well_depth, expected.intra_well_depth, 1e-8 * expected.intra_well_depth);
	EXPECT_NEAR(epithelial.inter_well_depth, 0.4 * epithelial.intra_well_depth, 1e-12 * expected.intra_well_depth);
	EXPECT_DOUBLE_EQ(epithelial.drag_coefficient, 5000.0 / 128.0);

	// no <motility>: the cells do not crawl
	EXPECT_FALSE(epithelial.is_motile);
}

TEST_F(ConfigReaderTest, ElementsGivenDirectlyAreTakenAsTheyAre)
{
	const simulation_config config = parse_simulation_config(write(settings {}));
	const cell_type& fibroblast = config.cell_types[1];

	EXPECT_EQ(fibroblast.name, "fibroblast");
	EXPECT_EQ(fibroblast.elements_per_cell, index_t { 19 });
	EXPECT_DOUBLE_EQ(fibroblast.equilibrium_distance, 1.5);
	EXPECT_DOUBLE_EQ(fibroblast.intra_well_depth, 2.0);
	EXPECT_DOUBLE_EQ(fibroblast.inter_well_depth, 0.6);
	EXPECT_DOUBLE_EQ(fibroblast.drag_coefficient, 3.0);

	EXPECT_TRUE(fibroblast.is_motile);
	EXPECT_DOUBLE_EQ(fibroblast.migration_speed, 0.7);
	EXPECT_DOUBLE_EQ(fibroblast.persistence_time, 5.0);
	EXPECT_DOUBLE_EQ(fibroblast.migration_bias, 0.25);
}

TEST_F(ConfigReaderTest, AffinitiesAreReadByNameAndAreOneWhereNotListed)
{
	const simulation_config config = parse_simulation_config(write(settings {}));

	// epithelial lists only the fibroblast (which is defined further down in the file)
	EXPECT_EQ(config.cell_types[0].adhesion_affinities, (std::vector<real_t> { 1.0, 0.3 }));
	EXPECT_EQ(config.cell_types[1].adhesion_affinities, (std::vector<real_t> { 0.2, 0.9 }));

	// the pair uses the geometric mean of the two entries
	const interaction_table table(config.cell_types, config.params);
	EXPECT_NEAR(table.inter(0, 1).attraction_scale, std::sqrt(0.3 * 0.2), 1e-15);
}

TEST_F(ConfigReaderTest, CellTypesAreOrderedByTheirIdNotByTheirPlaceInTheFile)
{
	settings s;
	replace(s.fibroblast, R"(ID="1")", R"(ID="0")");
	std::string text = s.text();
	replace(text, R"(name="epithelial" ID="0")", R"(name="epithelial" ID="1")");

	const simulation_config config = parse_simulation_config(write("settings.xml", text));

	EXPECT_EQ(config.cell_types[0].name, "fibroblast");
	EXPECT_EQ(config.cell_types[1].name, "epithelial");
	// affinities follow the IDs: the fibroblast's entry for itself is now the first
	EXPECT_EQ(config.cell_types[0].adhesion_affinities, (std::vector<real_t> { 0.9, 0.2 }));
}

TEST_F(ConfigReaderTest, OptionalPartsHaveDefaults)
{
	settings s;
	s.options.clear();
	s.micromechanics.clear();
	s.epithelial_affinities.clear();
	replace(s.epithelial_elements, "<adhesion_ratio>0.4</adhesion_ratio>", "");

	const simulation_config config = parse_simulation_config(write(s));

	EXPECT_FALSE(config.params.domain_walls);
	EXPECT_DOUBLE_EQ(config.params.steepness, 2.0);
	EXPECT_DOUBLE_EQ(config.params.cutoff_factor, 2.5);
	EXPECT_DOUBLE_EQ(config.params.thermal_energy, 0.0);
	EXPECT_EQ(config.cell_types[0].adhesion_affinities, (std::vector<real_t> { 1.0, 1.0 }));
	EXPECT_NEAR(config.cell_types[0].inter_well_depth, 0.5 * config.cell_types[0].intra_well_depth, 1e-12);
}

TEST_F(ConfigReaderTest, ReadsWhereTheCellPositionsAre)
{
	settings s;
	s.initial_conditions = R"(
    <initial_conditions>
        <cell_positions type="csv" enabled="true">
            <folder>./config</folder>
            <filename>cells.csv</filename>
        </cell_positions>
    </initial_conditions>)";
	EXPECT_EQ(parse_simulation_config(write(s)).cell_positions_file.value_or(std::filesystem::path {}),
			  std::filesystem::path("./config") / "cells.csv");

	replace(s.initial_conditions, "enabled=\"true\"", "enabled=\"false\"");
	EXPECT_FALSE(parse_simulation_config(write(s)).cell_positions_file.has_value());
}

// -------------------------------------------------------------------------------- what is refused

TEST_F(ConfigReaderTest, RefusesAFileThatIsMissingOrNotASettingsFile)
{
	EXPECT_THROW(parse_simulation_config(folder / "nowhere.xml"), std::runtime_error);
	EXPECT_THROW(parse_simulation_config(write("other.xml", "<?xml version=\"1.0\"?><something/>")),
				 std::runtime_error);
	EXPECT_THROW(parse_simulation_config(write("broken.xml", "<PhysiCell_settings>")), std::runtime_error);
}

TEST_F(ConfigReaderTest, RefusesMissingRequiredPartsAndSaysWhich)
{
	settings no_domain;
	no_domain.domain.clear();
	EXPECT_NE(error_of(no_domain).find("<domain>"), std::string::npos);

	settings no_timestep;
	replace(no_timestep.overall, "<dt_mechanics units=\"min\">0.0001</dt_mechanics>", "");
	EXPECT_NE(error_of(no_timestep).find("<dt_mechanics>"), std::string::npos);

	settings no_elements;
	no_elements.epithelial_elements.clear();
	const std::string message = error_of(no_elements);
	EXPECT_NE(message.find("epithelial"), std::string::npos) << message;
	EXPECT_NE(message.find("<sub_cellular_elements>"), std::string::npos) << message;

	settings no_stiffness;
	replace(no_stiffness.epithelial_elements, "<stiffness>5000</stiffness>", "");
	EXPECT_NE(error_of(no_stiffness).find("<stiffness>"), std::string::npos);
}

TEST_F(ConfigReaderTest, RefusesValuesThatMakeNoSenseAndNamesTheCellType)
{
	const auto with = [](void (*change)(settings&)) {
		settings s;
		change(s);
		return s;
	};

	EXPECT_NE(error_of(with([](settings& s) {
				  replace(s.epithelial_elements, "<elements>128", "<elements>0");
			  })).find("epithelial"),
			  std::string::npos);
	EXPECT_NE(error_of(with([](settings& s) {
				  replace(s.epithelial_elements, "<stiffness>5000", "<stiffness>-1");
			  })).find("epithelial"),
			  std::string::npos);
	EXPECT_NE(error_of(with([](settings& s) {
				  replace(s.fibroblast, "<drag_coefficient>3", "<drag_coefficient>0");
			  })).find("fibroblast"),
			  std::string::npos);
	EXPECT_NE(error_of(with([](settings& s) {
				  replace(s.fibroblast, "name=\"epithelial\">0.2", "name=\"epithelial\">1.2");
			  })).find("between 0 and 1"),
			  std::string::npos);
	EXPECT_NE(error_of(with([](settings& s) {
				  replace(s.fibroblast, "<persistence_time>5", "<persistence_time>0");
			  })).find("persistence_time"),
			  std::string::npos);
	EXPECT_NE(error_of(with([](settings& s) {
				  replace(s.micromechanics, "<steepness>3", "<steepness>0");
			  })).find("steepness"),
			  std::string::npos);
	EXPECT_NE(error_of(with([](settings& s) { replace(s.overall, ">0.0001<", ">0<"); })).find("dt_mechanics"),
			  std::string::npos);
	EXPECT_NE(error_of(with([](settings& s) { replace(s.domain, "<x_max>40", "<x_max>-50"); })).find("<domain>"),
			  std::string::npos);
}

TEST_F(ConfigReaderTest, RefusesAnAffinityToAnUnknownCellTypeAndBrokenIds)
{
	settings unknown;
	replace(unknown.epithelial_affinities, "name=\"fibroblast\"", "name=\"neuron\"");
	EXPECT_NE(error_of(unknown).find("neuron"), std::string::npos);

	settings gap;
	replace(gap.fibroblast, "ID=\"1\"", "ID=\"5\"");
	EXPECT_NE(error_of(gap).find("IDs"), std::string::npos);

	settings twice;
	replace(twice.fibroblast, "ID=\"1\"", "ID=\"0\"");
	EXPECT_NE(error_of(twice).find("IDs"), std::string::npos);

	settings same_name;
	replace(same_name.fibroblast, "name=\"fibroblast\" ID", "name=\"epithelial\" ID");
	EXPECT_NE(error_of(same_name).find("name"), std::string::npos);
}

TEST_F(ConfigReaderTest, InTwoDimensionsElementsMustBeGivenDirectly)
{
	settings flat;
	replace(flat.domain, "<use_2D>false</use_2D>", "<use_2D>true</use_2D>");

	// the epithelial type is given by stiffness and damping: that scaling is for three dimensions
	const std::string message = error_of(flat);
	EXPECT_NE(message.find("epithelial"), std::string::npos) << message;
	EXPECT_NE(message.find("<element>"), std::string::npos) << message;

	// with its elements given directly the same file is accepted
	flat.epithelial_elements = R"(
                    <sub_cellular_elements>
                        <elements>37</elements>
                        <element>
                            <equilibrium_distance>1</equilibrium_distance>
                            <intra_well_depth>1</intra_well_depth>
                            <inter_well_depth>0.5</inter_well_depth>
                            <drag_coefficient>1</drag_coefficient>
                        </element>
                    </sub_cellular_elements>)";
	const simulation_config config = parse_simulation_config(write(flat));
	EXPECT_EQ(config.dims, index_t { 2 });
	EXPECT_EQ(config.cell_types[0].elements_per_cell, index_t { 37 });
}

// ------------------------------------------------------------------------- from file to environment

TEST_F(ConfigReaderTest, EnvironmentFromAFileIsTheOneBuiltByHand)
{
	settings s;
	std::filesystem::create_directories(folder / "config");
	write("config/cells.csv", "x,y,z,type\n-12,0,0,0\n12,1,-2,fibroblast\n14.5,1,-2,1\n");
	s.initial_conditions = "\n    <initial_conditions>\n        <cell_positions type=\"csv\" enabled=\"true\">\n"
						   "            <folder>"
						   + (folder / "config").string() + "</folder>\n            <filename>cells.csv</filename>\n"
						   + "        </cell_positions>\n    </initial_conditions>";

	auto from_file = environment::create_from_config(write(s));

	// the same by hand
	const simulation_config config = parse_simulation_config(write(s));
	environment_builder builder;
	builder.set_time_step(0.0001);
	builder.resize(3, { -40, -30, -20 }, { 40, 30, 20 }, { 20, 20, 20 });
	builder.set_parameters(config.params);
	builder.add_cell_type(config.cell_types[0]);
	builder.add_cell_type(config.cell_types[1]);
	auto by_hand = builder.build();
	const std::array<real_t, 3> first = { -12.0, 0.0, 0.0 };
	const std::array<real_t, 3> second = { 12.0, 1.0, -2.0 };
	const std::array<real_t, 3> third = { 14.5, 1.0, -2.0 };
	by_hand->create_cell_of_type(0, first);
	by_hand->create_cell_of_type(1, second);
	by_hand->create_cell_of_type(1, third);
	by_hand->solver->initialize(*by_hand);

	ASSERT_EQ(from_file->cells->size(), 3U);
	EXPECT_EQ(from_file->agents->size(), 128U + 19U + 19U);
	EXPECT_EQ(from_file->cells->get_agent_at(0)->type_index(), index_t { 0 });
	EXPECT_EQ(from_file->cells->get_agent_at(1)->type_index(), index_t { 1 });
	EXPECT_EQ(from_file->cells->get_agent_at(2)->type_index(), index_t { 1 });
	EXPECT_EQ(from_file->cells->get_agent_at(1)->is_motile(), 1);
	EXPECT_NEAR(from_file->cells->get_agent_at(1)->position()[0], 12.0, 1e-9);

	// and they run alike, bit for bit
	for (auto* env : { from_file.get(), by_hand.get() })
	{
		env->solver->solve(*env, 50);
	}
	ASSERT_EQ(from_file->agents->size(), by_hand->agents->size());
	for (index_t i = 0; i < from_file->agents->size(); ++i)
	{
		for (std::size_t d = 0; d < 3; ++d)
		{
			ASSERT_EQ(from_file->agents->get_agent_at(i)->position()[d],
					  by_hand->agents->get_agent_at(i)->position()[d])
				<< "sub-agent " << i << " axis " << d;
		}
	}
	EXPECT_EQ(from_file->cells->get_agent_at(0)->pressure(), by_hand->cells->get_agent_at(0)->pressure());
}

TEST_F(ConfigReaderTest, RefusesCellPositionsThatCannotBeUsed)
{
	std::filesystem::create_directories(folder / "config");
	settings s;
	s.initial_conditions = "\n    <initial_conditions>\n        <cell_positions type=\"csv\" enabled=\"true\">\n"
						   "            <folder>"
						   + (folder / "config").string() + "</folder>\n            <filename>cells.csv</filename>\n"
						   + "        </cell_positions>\n    </initial_conditions>";
	const auto config_file = write(s);

	// no such file
	EXPECT_THROW(environment::create_from_config(config_file), std::runtime_error);

	const auto error_with = [&](const std::string& cells) {
		write("config/cells.csv", cells);
		try
		{
			environment::create_from_config(config_file);
		}
		catch (const std::runtime_error& e)
		{
			return std::string(e.what());
		}
		return std::string {};
	};

	EXPECT_NE(error_with("0,0,0,neuron\n").find("neuron"), std::string::npos);
	EXPECT_NE(error_with("0,0,0,7\n").find("unknown cell type"), std::string::npos);
	EXPECT_NE(error_with("0,0,0,0\n1,2\n").find("line 2"), std::string::npos);
	// the domain has walls and ends at x = 40
	EXPECT_NE(error_with("60,0,0,1\n").find("outside the domain"), std::string::npos);

	// a header, blank lines and further columns are fine
	EXPECT_TRUE(error_with("x,y,z,type,volume\n\n0,0,0,1,523.6\n").empty());
}

} // namespace physicore::mechanics::micromechanics::tests
