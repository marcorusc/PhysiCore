#pragma once

#include <span>
#include <vector>

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

	/**
	 * @brief Grows the cell at cell_index by one sub-agent.
	 *
	 * Several candidate positions are tried, each one equilibrium distance from a randomly chosen
	 * sub-agent of the cell in a random direction, and the one with the most room around it is used.
	 * The new sub-agent therefore never starts on top of another one. The random choices depend only
	 * on params.random_seed and on the number of growth events so far.
	 *
	 * Element parameters are not changed: a cell with twice the sub-agents has twice the volume.
	 */
	agent_interface* grow_cell(index_t cell_index);

	/**
	 * @brief Divides the cell at cell_index in two; no sub-agent moves.
	 *
	 * A plane perpendicular to normal splits the sub-agents. The daughter_fraction of them that lies
	 * furthest along normal becomes a new cell, which is returned; the rest stays in the original cell.
	 * 0.5 divides symmetrically; 0.3 gives daughters holding 30% and 70% of the sub-agents. Each
	 * daughter keeps at least one sub-agent. The new cell copies the properties of the original.
	 *
	 * @throws std::invalid_argument if normal has the wrong size or is zero, daughter_fraction is not
	 *         between 0 and 1, or the cell has fewer than two sub-agents
	 */
	cell_interface* divide_cell(index_t cell_index, std::span<const real_t> normal, real_t daughter_fraction = 0.5);

	/**
	 * @brief Unit vector along which the cell at cell_index is most stretched.
	 *
	 * An axis has no direction of its own; the vector is returned with its largest component positive.
	 * For a cell with fewer than two sub-agents, or one that is not stretched at all, the result is a
	 * valid unit vector without further meaning.
	 */
	std::vector<real_t> longest_axis(index_t cell_index);

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

private:
	// number of grow_cell() calls so far; with the seed it identifies the random choices of a call
	index_t growth_events = 0;
};

} // namespace physicore::mechanics::micromechanics
