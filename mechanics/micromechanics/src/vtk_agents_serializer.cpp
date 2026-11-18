#include "vtk_agents_serializer.h"

#include <array>
#include <filesystem>
#include <iomanip>
#include <memory>
#include <sstream>
#include <vtkCellArray.h>
#include <vtkCellType.h>
#include <vtkPointData.h>
#include <vtkPoints.h>

#include "micromechanics/vtk_serializer.h"

namespace physicore::mechanics::micromechanics {

namespace {
// Value of sub-agent i padded to three components, as VTK expects for points and vectors
std::array<double, 3> padded(const std::vector<real_t>& values, index_t i, index_t dims)
{
	std::array<double, 3> result = { 0.0, 0.0, 0.0 };
	for (index_t d = 0; d < dims && d < 3; ++d)
	{
		result.at(d) = values[i * dims + d];
	}
	return result;
}
} // namespace

serializer_ptr make_vtk_agents_serializer(std::string_view output_dir)
{
	return std::make_unique<vtk_agents_serializer>(output_dir);
}

vtk_agents_serializer::vtk_agents_serializer(std::string_view output_dir)
	: vtk_serializer_base(output_dir, "vtk_micromechanics_agents", "micromechanics_agents.pvd")
{
	cell_index_array->SetName("cell_index");
	cell_index_array->SetNumberOfComponents(1);
	unstructured_grid->GetPointData()->AddArray(cell_index_array);

	cell_type_array->SetName("cell_type");
	cell_type_array->SetNumberOfComponents(1);
	unstructured_grid->GetPointData()->AddArray(cell_type_array);

	velocity_array->SetName("velocity");
	velocity_array->SetNumberOfComponents(3);
	unstructured_grid->GetPointData()->AddArray(velocity_array);

	writer->SetInputData(unstructured_grid);
	writer->SetCompressorTypeToNone();
}

void vtk_agents_serializer::serialize(const environment& e, real_t current_time)
{
	const agent_data& data = generic_agent_solver<agent>::retrieve_agent_data(*e.agents);
	const cell_data& cells = generic_agent_solver<cell>::retrieve_agent_data(*e.cells);

	const index_t dims = data.base_data.dims;
	const auto agents_count = static_cast<vtkIdType>(data.agents_count);

	auto points = vtkSmartPointer<vtkPoints>::New();
	points->SetNumberOfPoints(agents_count);

	auto vertices = vtkSmartPointer<vtkCellArray>::New();

	cell_index_array->SetNumberOfTuples(agents_count);
	cell_type_array->SetNumberOfTuples(agents_count);
	velocity_array->SetNumberOfTuples(agents_count);

	for (vtkIdType id = 0; id < agents_count; ++id)
	{
		const auto i = static_cast<index_t>(id);
		const index_t cell_index = data.cell_indices[i];
		const bool has_cell = cell_index < cells.cells_count;

		points->SetPoint(id, padded(data.base_data.positions, i, dims).data());

		vertices->InsertNextCell(1);
		vertices->InsertCellPoint(id);

		cell_index_array->SetValue(id, has_cell ? static_cast<vtkIdType>(cell_index) : -1);
		cell_type_array->SetValue(id, has_cell ? static_cast<vtkIdType>(cells.type_indices[cell_index]) : -1);
		velocity_array->SetTuple(id, padded(data.velocities, i, dims).data());
	}

	unstructured_grid->SetPoints(points);
	unstructured_grid->SetCells(VTK_VERTEX, vertices);

	std::ostringstream ss;
	ss << "micromechanics_agents_" << std::setw(6) << std::setfill('0') << iteration << ".vtu";

	const auto file_name = ss.str();
	const auto file_path = std::filesystem::path(vtks_dir) / file_name;

	writer->SetFileName(file_path.string().c_str());
	writer->Write();

	append_to_pvd(file_name, current_time);

	iteration++;
}

} // namespace physicore::mechanics::micromechanics
