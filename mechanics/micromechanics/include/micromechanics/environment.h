#pragma once

#include <span>
#include <vector>

#include <common/generic_agent_solver.h>
#include <common/mesh.h>
#include <common/timestep_executor.h>
#include <common/types.h>

#include "agent_container.h"
#include "cell_container.h"
#include "cell_type.h"
#include "serializer.h"
#include "simulation_parameters.h"
#include "solver.h"

namespace physicore::mechanics::micromechanics {

/**
 * @brief Sub-cellular element model environment: cells, their sub-agents and the solver acting on them.
 *
 * Every sub-agent refers to its cell by the cell's index in the cell container. Containers remove
 * elements by moving the last one into the freed slot, so cells and sub-agents must be created and
 * removed through this class, which keeps the agent -> cell indices and the per-cell agent counts valid.
 */
class environment : public timestep_executor,
					private physicore::generic_agent_solver<agent>,
					private physicore::generic_agent_solver<cell>
{
public:
	environment(const cartesian_mesh& mesh, real_t timestep);

	void run_single_timestep() override;

	void serialize_state(real_t current_time) override;

	// Adds a cell type and returns its index. There is always the type 0, which starts with the
	// default values of cell_type and can be changed in place (cell_types[0]).
	index_t add_cell_type(cell_type type);

	// Creates a cell of the given type without sub-agents; its index is cells->size() - 1 until a
	// cell is removed.
	// @throws std::out_of_range if there is no such type
	cell_interface* create_cell(index_t type_index = 0);

	// Creates a cell of type 0 of elements_count sub-agents packed around centre without overlaps,
	// one equilibrium distance of the type apart (see packed_ball_positions). centre holds one value
	// per dimension.
	cell_interface* create_cell(std::span<const real_t> centre, index_t elements_count);

	// The same for a cell of the given type, with the number of sub-agents of that type
	// (cell_type::elements_per_cell)
	cell_interface* create_cell_of_type(index_t type_index, std::span<const real_t> centre);

	// The same with a number of sub-agents of one's own choosing
	cell_interface* create_cell_of_type(index_t type_index, std::span<const real_t> centre, index_t elements_count);

	// Creates a sub-agent that belongs to the cell at cell_index
	agent_interface* create_agent(index_t cell_index);

	// Removes a single sub-agent; the last sub-agent takes over its index
	void remove_agent(index_t agent_index);

	// Removes a cell together with all its sub-agents; the last cell takes over its index
	void remove_cell(index_t cell_index);

	cell_container_ptr cells;
	agent_container_ptr agents;
	solver_ptr solver;
	serializer_ptr serializer;

	real_t mechanics_timestep;
	cartesian_mesh mesh;
	simulation_parameters params;

	// The kinds of cells; a cell refers to one by its index (cell_interface::type_index)
	std::vector<cell_type> cell_types = { cell_type {} };
};

} // namespace physicore::mechanics::micromechanics
