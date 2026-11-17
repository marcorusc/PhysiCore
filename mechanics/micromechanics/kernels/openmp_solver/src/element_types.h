#pragma once

#include <common/types.h>

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

// Cell type of a sub-agent: that of its cell. A sub-agent without a cell is treated as type 0.
inline index_t type_of_element(index_t cell_index, const index_t* type_indices, index_t cells_count)
{
	return cell_index < cells_count ? type_indices[cell_index] : 0;
}

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
