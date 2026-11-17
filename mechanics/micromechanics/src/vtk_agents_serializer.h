#pragma once

#include <string_view>
#include <vtkIdTypeArray.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <common/generic_agent_solver.h>
#include <common/types.h>
#include <common/vtk_serializer_base.h>

#include "micromechanics/environment.h"
#include "micromechanics/serializer.h"

namespace physicore::mechanics::micromechanics {

class vtk_agents_serializer : public vtk_serializer_base,
							  public serializer,
							  private generic_agent_solver<agent>,
							  private generic_agent_solver<cell>
{
	vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();

	vtkSmartPointer<vtkIdTypeArray> cell_index_array = vtkSmartPointer<vtkIdTypeArray>::New();
	vtkSmartPointer<vtkRealArray> agent_radius_array = vtkSmartPointer<vtkRealArray>::New();
	vtkSmartPointer<vtkRealArray> velocity_array = vtkSmartPointer<vtkRealArray>::New();
	vtkSmartPointer<vtkRealArray> previous_velocity_array = vtkSmartPointer<vtkRealArray>::New();

public:
	explicit vtk_agents_serializer(std::string_view output_dir);

	void serialize(const environment& e, real_t current_time) override;
};

} // namespace physicore::mechanics::micromechanics
