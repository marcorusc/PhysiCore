#pragma once

#include <memory>

#include <common/types.h>

namespace physicore::mechanics::micromechanics {

class environment;

/**
 * @brief Contract every micromechanics backend (kernel) implements.
 *
 * Backends own all hardware-specific work: neighbor search, force evaluation,
 * per-cell aggregation and position integration.
 */
class solver
{
public:
	solver() = default;
	solver(const solver&) = delete;
	solver& operator=(const solver&) = delete;
	solver(solver&&) = delete;
	solver& operator=(solver&&) = delete;

	// Prepare backend state for the given environment
	virtual void initialize(environment& e) = 0;

	// Advance the mechanics by the given number of mechanics timesteps
	virtual void solve(environment& e, index_t iterations) = 0;

	virtual ~solver() = default;
};

using solver_ptr = std::unique_ptr<solver>;

} // namespace physicore::mechanics::micromechanics
