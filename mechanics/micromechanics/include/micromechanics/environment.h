#pragma once

#include <span>

#include <common/generic_agent_solver.h>
#include <common/mesh.h>
#include <common/timestep_executor.h>
#include <common/types.h>

#include "agent_container.h"
#include "cell_container.h"
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

	// Creates a cell without sub-agents; its index is cells->size() - 1 until a cell is removed.
	// Its sub-agent radius is half the current equilibrium distance, so sub-agents drawn as
	// spheres touch at equilibrium. Set params before creating cells.
	cell_interface* create_cell();

	// Creates a cell of elements_count sub-agents packed around centre without overlaps, one
	// equilibrium distance apart (see packed_ball_positions). centre holds one value per dimension.
	cell_interface* create_cell(std::span<const real_t> centre, index_t elements_count);

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
};

} // namespace physicore::mechanics::micromechanics
