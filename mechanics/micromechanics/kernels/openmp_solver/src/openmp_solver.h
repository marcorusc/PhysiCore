#pragma once

#include <micromechanics/solver.h>

#include "force_solver.h"
#include "position_solver.h"

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

class openmp_solver : public solver
{
	force_solver forces;
	position_solver positions;

public:
	void initialize(environment& e) override;
	void solve(environment& e, index_t iterations) override;
};

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
