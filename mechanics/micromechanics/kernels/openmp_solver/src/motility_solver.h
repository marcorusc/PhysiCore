#pragma once

#include <common/generic_agent_solver.h>
#include <common/types.h>
#include <micromechanics/environment.h>

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

/**
 * @brief Chooses the direction in which every motile cell crawls: a persistent random walk.
 *
 * A motile cell keeps its motility direction for a while and then turns. In every step it turns with
 * probability timestep / persistence_time, so the straight runs last persistence_time on average. The
 * new direction is
 *
 *     normalise((1 - migration_bias) * xi + migration_bias * migration_bias_direction)
 *
 * with xi a random unit vector, uniform on the circle or sphere. A bias of 0 gives an unbiased walk, a
 * bias of 1 a cell that always follows its bias direction. This is the motility model of PhysiCell
 * (Ghaffarizadeh et al. 2018), with the same parameters.
 *
 * A motile cell without a direction (all zeros, as after creation) draws one in its first step. If the
 * mix has no length, the cell keeps the direction it had.
 *
 * The solver only sets motility_direction. position_solver turns it into motion by adding
 * migration_speed * motility_direction to the velocity of every sub-agent of the cell.
 *
 * The random numbers are a pure function of (random_seed, step, cell index), so a run does not depend
 * on the number of threads, and they are independent of the thermal noise.
 */
class motility_solver : private generic_agent_solver<cell>
{
	// number of updates done so far; together with the seed it identifies the random numbers of a step
	index_t step_ = 0;

public:
	// Throws std::invalid_argument if a motile cell has a persistence time that is not positive or a
	// migration bias outside [0, 1]
	void validate(environment& e);

	void update(environment& e);
};

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
