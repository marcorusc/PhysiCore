#pragma once

#include <cmath>
#include <cstdint>
#include <numbers>

#include <common/types.h>

namespace physicore::mechanics::micromechanics::kernels::openmp_solver {

/**
 * Counter-based random numbers: a value is a pure function of (seed, step, sub-agent, component) and
 * not the next output of a shared generator. A simulation is therefore reproducible bit for bit from
 * its seed, whatever the number of threads or the order in which sub-agents are processed.
 */

// Finalizer of the SplitMix64 generator: a bijection on 64-bit integers that mixes all input bits
constexpr std::uint64_t mix(std::uint64_t z) noexcept
{
	z += 0x9E3779B97F4A7C15ULL;
	z = (z ^ (z >> 30U)) * 0xBF58476D1CE4E5B9ULL;
	z = (z ^ (z >> 27U)) * 0x94D049BB133111EBULL;
	return z ^ (z >> 31U);
}

// 64 random bits for the given counters; stream tells apart several values needed for the same counters
constexpr std::uint64_t random_bits(std::uint64_t seed, std::uint64_t step, std::uint64_t agent,
									std::uint64_t component, std::uint64_t stream) noexcept
{
	std::uint64_t h = mix(seed);
	h = mix(h ^ step);
	h = mix(h ^ agent);
	h = mix(h ^ (2 * component + stream));
	return h;
}

// Uniformly distributed in (0, 1]: never zero, so its logarithm is finite
constexpr real_t to_unit_interval(std::uint64_t bits) noexcept
{
	return (static_cast<real_t>(bits >> 11U) + 1.0) * 0x1.0p-53;
}

// Normally distributed with zero mean and unit variance (Box-Muller transform)
inline real_t standard_normal(std::uint64_t seed, std::uint64_t step, std::uint64_t agent,
							  std::uint64_t component) noexcept
{
	const real_t u1 = to_unit_interval(random_bits(seed, step, agent, component, 0));
	const real_t u2 = to_unit_interval(random_bits(seed, step, agent, component, 1));

	return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * std::numbers::pi * u2);
}

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver
