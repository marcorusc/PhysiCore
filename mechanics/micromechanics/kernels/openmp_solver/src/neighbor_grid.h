#pragma once

#include <array>
#include <vector>

#include <common/types.h>

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

/**
 * @brief Uniform grid that sorts sub-agents into bins at least as wide as the interaction range.
 *
 * Two sub-agents closer than the interaction range are then always in the same or in adjacent bins,
 * so a sub-agent only has to look at 3^dims bins instead of at every other sub-agent.
 *
 * The grid covers the bounding box of the sub-agents and is rebuilt from scratch from the current
 * positions. Its memory is bounded: when the sub-agents are spread too thinly the bins are made wider,
 * which stays correct and only costs speed.
 */
class neighbor_grid
{
	index_t dims_ = 0;
	real_t inverse_bin_size_ = 0.0;
	std::array<real_t, 3> origin_ {};
	std::array<index_t, 3> shape_ { 1, 1, 1 };

	// agents of bin b are sorted_agents_[bin_starts_[b]] ... sorted_agents_[bin_starts_[b + 1] - 1]
	std::vector<index_t> bin_starts_;
	std::vector<index_t> sorted_agents_;
	std::vector<index_t> bin_of_agent_;

	std::array<index_t, 3> bin_coordinates(const real_t* position) const;

public:
	// positions holds agents_count * dims values; interaction_range must be positive
	// @throws std::runtime_error if a position is not a finite number
	void build(const real_t* positions, index_t agents_count, index_t dims, real_t interaction_range);

	// Number of bins along each axis (1 for axes beyond dims)
	const std::array<index_t, 3>& shape() const { return shape_; }

	/**
	 * @brief Calls visit(j) for every sub-agent j in the bin containing position and in the bins around it.
	 *
	 * The sub-agent at position itself is among them, and so are sub-agents further away than the
	 * interaction range: the caller filters by distance. The order of the calls depends only on the
	 * positions the grid was built from.
	 */
	template <typename Visitor>
	void visit_candidates(const real_t* position, const Visitor& visit) const
	{
		const std::array<index_t, 3> centre = bin_coordinates(position);

		std::array<index_t, 3> first {};
		std::array<index_t, 3> last {};
		for (std::size_t d = 0; d < 3; ++d)
		{
			first[d] = centre[d] == 0 ? 0 : centre[d] - 1;
			last[d] = centre[d] + 1 < shape_[d] ? centre[d] + 1 : shape_[d] - 1;
		}

		for (index_t z = first[2]; z <= last[2]; ++z)
		{
			for (index_t y = first[1]; y <= last[1]; ++y)
			{
				for (index_t x = first[0]; x <= last[0]; ++x)
				{
					const index_t bin = x + shape_[0] * (y + shape_[1] * z);

					for (index_t k = bin_starts_[bin]; k < bin_starts_[bin + 1]; ++k)
					{
						visit(sorted_agents_[k]);
					}
				}
			}
		}
	}
};

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
