#include "micromechanics/cell_type.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace physicore::mechanics::micromechanics {

namespace {
void check(const cell_type& type, const simulation_parameters& params)
{
	if (!(params.steepness > 0.0) || !(params.cutoff_factor > 0.0))
	{
		throw std::invalid_argument("Steepness and cutoff factor must be positive");
	}
	if (!(type.drag_coefficient > 0.0) || !(type.equilibrium_distance > 0.0) || !(type.intra_well_depth > 0.0))
	{
		throw std::invalid_argument("Cell type '" + type.name
									+ "': drag coefficient, equilibrium distance and intra-cellular well depth "
									  "must be positive");
	}
	if (!(type.inter_well_depth >= 0.0))
	{
		throw std::invalid_argument("Cell type '" + type.name
									+ "': the inter-cellular well depth must not be negative");
	}
	for (const real_t affinity : type.adhesion_affinities)
	{
		if (std::isnan(affinity) || affinity < 0.0 || affinity > 1.0)
		{
			throw std::invalid_argument("Cell type '" + type.name + "': adhesion affinities must be between 0 and 1");
		}
	}
}
} // namespace

interaction_table::interaction_table(std::span<const cell_type> types, const simulation_parameters& params)
	: types_count_(types.size())
{
	if (types.empty())
	{
		throw std::invalid_argument("At least one cell type is needed");
	}

	intra_.reserve(types_count_);
	inter_.reserve(types_count_ * types_count_);

	for (index_t t = 0; t < types_count_; ++t)
	{
		check(types[t], params);
		intra_.push_back({ .potential = intra_potential(types[t], params), .attraction_scale = 1.0 });
		largest_cutoff_ = std::max(largest_cutoff_, intra_.back().potential.cutoff_distance);
	}

	for (index_t t = 0; t < types_count_; ++t)
	{
		for (index_t s = 0; s < types_count_; ++s)
		{
			if (t == s)
			{
				// the general rule below gives exactly this; written out so that no rounding of the
				// square roots enters when there is one type only
				inter_.push_back(
					{ .potential = inter_potential(types[t], params), .attraction_scale = types[t].affinity_to(t) });
				continue;
			}

			const real_t equilibrium_distance = 0.5 * (types[t].equilibrium_distance + types[s].equilibrium_distance);
			// both products are symmetric in t and s, so the entries (t, s) and (s, t) are identical
			const real_t well_depth = std::sqrt(types[t].inter_well_depth * types[s].inter_well_depth);
			const real_t affinity = std::sqrt(types[t].affinity_to(s) * types[s].affinity_to(t));

			inter_.push_back({ .potential = { .well_depth = well_depth,
											  .equilibrium_distance = equilibrium_distance,
											  .steepness = params.steepness,
											  .cutoff_distance = params.cutoff_factor * equilibrium_distance },
							   .attraction_scale = affinity });
			largest_cutoff_ = std::max(largest_cutoff_, inter_.back().potential.cutoff_distance);
		}
	}
}

} // namespace physicore::mechanics::micromechanics
