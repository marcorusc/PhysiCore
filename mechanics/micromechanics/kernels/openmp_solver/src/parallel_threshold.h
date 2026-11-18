#pragma once

#include <common/types.h>

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

// Loops over fewer items than this run on the calling thread: for small systems starting a parallel
// region costs more than the loop itself.
inline constexpr index_t parallel_threshold = 256;

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
