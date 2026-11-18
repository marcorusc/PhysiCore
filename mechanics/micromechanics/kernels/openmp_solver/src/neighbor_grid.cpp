#include "neighbor_grid.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

namespace {
// Upper bound on the number of bins along one axis, chosen so that the whole grid stays within a few
// bins per sub-agent plus a constant.
index_t max_bins_per_axis(index_t agents_count, index_t dims)
{
	const real_t max_bins = 4.0 * static_cast<real_t>(agents_count) + 4096.0;
	const real_t per_axis = std::floor(std::pow(max_bins, 1.0 / static_cast<real_t>(dims)));

	return std::max<index_t>(1, static_cast<index_t>(per_axis));
}
} // namespace

std::array<index_t, 3> neighbor_grid::bin_coordinates(const real_t* position) const
{
	std::array<index_t, 3> coordinates {};
	for (index_t d = 0; d < dims_; ++d)
	{
		// positions outside the box the grid was built for land in the nearest boundary bin
		const real_t scaled = std::max(0.0, (position[d] - origin_[d]) * inverse_bin_size_);
		coordinates[d] = std::min(shape_[d] - 1, static_cast<index_t>(scaled));
	}
	return coordinates;
}

void neighbor_grid::build(const real_t* positions, index_t agents_count, index_t dims, real_t interaction_range)
{
	assert(dims >= 1 && dims <= 3);
	assert(interaction_range > 0.0);

	dims_ = dims;
	origin_ = { 0.0, 0.0, 0.0 };
	shape_ = { 1, 1, 1 };

	std::array<real_t, 3> highest {};
	if (agents_count > 0)
	{
		for (index_t d = 0; d < dims; ++d)
		{
			origin_[d] = std::numeric_limits<real_t>::max();
			highest[d] = std::numeric_limits<real_t>::lowest();
		}
		bool finite = true;
		for (index_t i = 0; i < agents_count; ++i)
		{
			for (index_t d = 0; d < dims; ++d)
			{
				finite = finite && std::isfinite(positions[i * dims + d]);
				origin_[d] = std::min(origin_[d], positions[i * dims + d]);
				highest[d] = std::max(highest[d], positions[i * dims + d]);
			}
		}
		if (!finite)
		{
			// no grid can hold a sub-agent at infinity, and nothing that follows would mean anything
			throw std::runtime_error("A sub-agent position is no longer a finite number. The simulation has "
									 "blown up, most likely because the timestep is too large for the stiffness of "
									 "the cells (see largest_stable_timestep).");
		}
	}

	// bins must not be narrower than the interaction range; they get wider if there would be too many
	// the small margin keeps a pair at exactly the interaction range in adjacent bins despite rounding
	real_t bin_size = interaction_range * (1.0 + 1e-9);
	const index_t bins_limit = max_bins_per_axis(agents_count, dims);
	for (index_t d = 0; d < dims; ++d)
	{
		bin_size = std::max(bin_size, (highest[d] - origin_[d]) / static_cast<real_t>(bins_limit));
	}
	inverse_bin_size_ = 1.0 / bin_size;

	index_t bins_count = 1;
	for (index_t d = 0; d < dims; ++d)
	{
		const auto needed = static_cast<index_t>((highest[d] - origin_[d]) * inverse_bin_size_) + 1;
		shape_[d] = std::min(needed, bins_limit);
		bins_count *= shape_[d];
	}

	// counting sort of the sub-agents by bin
	bin_of_agent_.resize(agents_count);
	bin_starts_.assign(bins_count + 1, 0);
	for (index_t i = 0; i < agents_count; ++i)
	{
		const std::array<index_t, 3> c = bin_coordinates(positions + i * dims);
		bin_of_agent_[i] = c[0] + shape_[0] * (c[1] + shape_[1] * c[2]);
		++bin_starts_[bin_of_agent_[i] + 1];
	}
	for (index_t bin = 0; bin < bins_count; ++bin)
	{
		bin_starts_[bin + 1] += bin_starts_[bin];
	}

	sorted_agents_.resize(agents_count);
	std::vector<index_t> next_slot(bin_starts_.begin(), bin_starts_.end() - 1);
	for (index_t i = 0; i < agents_count; ++i)
	{
		sorted_agents_[next_slot[bin_of_agent_[i]]++] = i;
	}
}

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
