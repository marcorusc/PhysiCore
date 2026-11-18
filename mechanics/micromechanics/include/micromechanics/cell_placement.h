#pragma once

#include <vector>

#include <common/types.h>
#include <micromechanics/micromechanics_export.h>

namespace physicore::mechanics::micromechanics {

/**
 * @brief Positions of count elements packed as closely as possible around the origin.
 *
 * The elements sit on the densest regular lattice with nearest-neighbour distance spacing: a line in
 * one dimension, the hexagonal lattice in two and the face-centred cubic lattice in three. Of all
 * lattice sites the count closest to the origin are used, so the result is as round as a lattice
 * allows, no two elements are closer than spacing, and none overlaps.
 *
 * The positions are ordered by distance from the origin, starting with the origin itself, and the
 * positions for a smaller count are the beginning of those for a larger one.
 *
 * @return count * dims coordinates
 * @throws std::invalid_argument if dims is not 1, 2 or 3, or spacing is not positive
 */
MICROMECHANICS_EXPORT std::vector<real_t> packed_ball_positions(index_t dims, index_t count, real_t spacing);

} // namespace physicore::mechanics::micromechanics
