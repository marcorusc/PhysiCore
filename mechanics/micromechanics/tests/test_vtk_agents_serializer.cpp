#include <algorithm>
#include <array>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vtkIdTypeArray.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>

#include <common/vtk_serializer_base.h>
#include <gtest/gtest.h>

#include "micromechanics/environment.h"
#include "micromechanics/vtk_serializer.h"

namespace physicore::mechanics::micromechanics::tests {

namespace {
environment make_environment(index_t dims)
{
	return { { dims, { -100, -100, -100 }, { 100, 100, 100 }, { 20, 20, 20 } }, 0.1 };
}

vtkSmartPointer<vtkUnstructuredGrid> read_grid(const std::filesystem::path& file)
{
	auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
	reader->SetFileName(file.string().c_str());
	reader->Update();
	return reader->GetOutput();
}

std::array<double, 3> tuple_of(vtkPointData* point_data, const char* name, vtkIdType id)
{
	std::array<double, 3> tuple {};
	auto* array = vtkRealArray::SafeDownCast(point_data->GetArray(name));
	if (array == nullptr)
	{
		ADD_FAILURE() << "Missing array: " << name;
		return tuple;
	}
	array->GetTuple(id, tuple.data());
	return tuple;
}

void expect_tuple(const std::array<double, 3>& actual, const std::array<double, 3>& expected)
{
	for (std::size_t d = 0; d < 3; ++d)
	{
		EXPECT_DOUBLE_EQ(actual.at(d), expected.at(d)) << "component " << d;
	}
}

std::string contents_of(const std::filesystem::path& file)
{
	const std::ifstream stream(file);
	std::ostringstream contents;
	contents << stream.rdbuf();
	return contents.str();
}
} // namespace

class VtkMicromechanicsAgentsSerializerTest : public ::testing::Test
{
protected:
	std::filesystem::path output_dir;

	void SetUp() override
	{
		const auto* test = ::testing::UnitTest::GetInstance()->current_test_info();
		output_dir =
			std::filesystem::temp_directory_path() / (std::string("physicore_micromechanics_vtk_") + test->name());
		std::filesystem::remove_all(output_dir);
	}

	void TearDown() override { std::filesystem::remove_all(output_dir); }

	std::filesystem::path vtu(index_t iteration) const
	{
		std::ostringstream name;
		name << "micromechanics_agents_00000" << iteration << ".vtu";
		return output_dir / "vtk_micromechanics_agents" / name.str();
	}
};

TEST_F(VtkMicromechanicsAgentsSerializerTest, WritesPositionsAndAllArrays)
{
	environment env = make_environment(3);
	env.create_cell();
	env.create_cell(env.add_cell_type(cell_type {}));

	auto* first = env.create_agent(1);
	std::ranges::copy(std::array { 1.0, 2.0, 3.0 }, first->position().begin());
	std::ranges::copy(std::array { 0.1, 0.2, 0.3 }, first->velocity().begin());

	auto* second = env.create_agent(0);
	std::ranges::copy(std::array { -4.0, 5.0, -6.0 }, second->position().begin());

	// a sub-agent that belongs to no cell
	env.agents->create()->position()[0] = 9.0;

	env.serializer = make_vtk_agents_serializer(output_dir.string());
	env.serialize_state(0.0);

	ASSERT_TRUE(std::filesystem::exists(vtu(0)));
	const auto grid = read_grid(vtu(0));
	ASSERT_EQ(grid->GetNumberOfPoints(), 3);
	ASSERT_EQ(grid->GetNumberOfCells(), 3);

	std::array<double, 3> point {};
	grid->GetPoint(0, point.data());
	expect_tuple(point, { 1.0, 2.0, 3.0 });
	grid->GetPoint(1, point.data());
	expect_tuple(point, { -4.0, 5.0, -6.0 });
	grid->GetPoint(2, point.data());
	expect_tuple(point, { 9.0, 0.0, 0.0 });

	auto* point_data = grid->GetPointData();
	ASSERT_EQ(point_data->GetNumberOfArrays(), 3);

	auto* cell_index = vtkIdTypeArray::SafeDownCast(point_data->GetArray("cell_index"));
	ASSERT_NE(cell_index, nullptr);
	EXPECT_EQ(cell_index->GetValue(0), 1);
	EXPECT_EQ(cell_index->GetValue(1), 0);
	EXPECT_EQ(cell_index->GetValue(2), -1);

	auto* cell_type_values = vtkIdTypeArray::SafeDownCast(point_data->GetArray("cell_type"));
	ASSERT_NE(cell_type_values, nullptr);
	EXPECT_EQ(cell_type_values->GetValue(0), 1); // the second cell is of the second type
	EXPECT_EQ(cell_type_values->GetValue(1), 0);
	EXPECT_EQ(cell_type_values->GetValue(2), -1);

	expect_tuple(tuple_of(point_data, "velocity", 0), { 0.1, 0.2, 0.3 });
	expect_tuple(tuple_of(point_data, "velocity", 1), { 0.0, 0.0, 0.0 });
}

TEST_F(VtkMicromechanicsAgentsSerializerTest, TwoDimensionalDataIsPaddedWithZeros)
{
	environment env = make_environment(2);
	env.create_cell();
	auto* element = env.create_agent(0);
	std::ranges::copy(std::array { 1.5, -2.5 }, element->position().begin());
	std::ranges::copy(std::array { 0.25, 0.5 }, element->velocity().begin());

	env.serializer = make_vtk_agents_serializer(output_dir.string());
	env.serialize_state(0.0);

	const auto grid = read_grid(vtu(0));
	ASSERT_EQ(grid->GetNumberOfPoints(), 1);

	std::array<double, 3> point {};
	grid->GetPoint(0, point.data());
	expect_tuple(point, { 1.5, -2.5, 0.0 });
	expect_tuple(tuple_of(grid->GetPointData(), "velocity", 0), { 0.25, 0.5, 0.0 });
}

TEST_F(VtkMicromechanicsAgentsSerializerTest, EverySerializationAddsAFileAndACollectionEntry)
{
	environment env = make_environment(3);
	env.create_cell();
	auto* element = env.create_agent(0);

	env.serializer = make_vtk_agents_serializer(output_dir.string());
	env.serialize_state(0.0);
	element->position()[0] = 2.0;
	env.serialize_state(1.5);

	ASSERT_TRUE(std::filesystem::exists(vtu(0)));
	ASSERT_TRUE(std::filesystem::exists(vtu(1)));

	// the second file holds the state at the second call
	std::array<double, 3> point {};
	read_grid(vtu(1))->GetPoint(0, point.data());
	EXPECT_DOUBLE_EQ(point[0], 2.0);

	const std::string collection = contents_of(output_dir / "micromechanics_agents.pvd");
	EXPECT_NE(collection.find("micromechanics_agents_000000.vtu"), std::string::npos);
	EXPECT_NE(collection.find("micromechanics_agents_000001.vtu"), std::string::npos);
	EXPECT_NE(collection.find("timestep=\"1.5"), std::string::npos);
	EXPECT_NE(collection.find("</VTKFile>"), std::string::npos);
}

TEST_F(VtkMicromechanicsAgentsSerializerTest, EmptySystemWritesAnEmptyFile)
{
	environment env = make_environment(3);

	env.serializer = make_vtk_agents_serializer(output_dir.string());
	env.serialize_state(0.0);

	ASSERT_TRUE(std::filesystem::exists(vtu(0)));
	EXPECT_EQ(read_grid(vtu(0))->GetNumberOfPoints(), 0);
}

} // namespace physicore::mechanics::micromechanics::tests
