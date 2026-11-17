#include <cmath>
#include <cstdint>
#include <set>
#include <vector>

#include <gtest/gtest.h>

#include "counter_random.h"

namespace physicore::mechanics::micromechanics::kernels::openmp_solver::tests {

namespace {
struct moments
{
	real_t mean = 0.0;
	real_t variance = 0.0;
	real_t skewness = 0.0;
	real_t kurtosis = 0.0;
};

moments moments_of(const std::vector<real_t>& samples)
{
	const auto count = static_cast<real_t>(samples.size());

	moments result;
	for (const real_t x : samples)
	{
		result.mean += x / count;
	}
	for (const real_t x : samples)
	{
		const real_t centred = x - result.mean;
		result.variance += centred * centred / count;
		result.skewness += centred * centred * centred / count;
		result.kurtosis += centred * centred * centred * centred / count;
	}
	result.skewness /= std::pow(result.variance, 1.5);
	result.kurtosis /= result.variance * result.variance;
	return result;
}

real_t correlation(const std::vector<real_t>& a, const std::vector<real_t>& b)
{
	const moments ma = moments_of(a);
	const moments mb = moments_of(b);

	real_t covariance = 0.0;
	for (std::size_t i = 0; i < a.size(); ++i)
	{
		covariance += (a[i] - ma.mean) * (b[i] - mb.mean) / static_cast<real_t>(a.size());
	}
	return covariance / std::sqrt(ma.variance * mb.variance);
}

constexpr index_t samples_count = 200000;
// with 200000 samples the standard error of the mean and of a correlation is about 0.0022
constexpr real_t five_sigma = 0.0112;
} // namespace

TEST(CounterRandomTest, SameCountersGiveTheSameNumber)
{
	EXPECT_EQ(standard_normal(7, 3, 11, 2), standard_normal(7, 3, 11, 2));
	EXPECT_EQ(random_bits(7, 3, 11, 2, 1), random_bits(7, 3, 11, 2, 1));
}

TEST(CounterRandomTest, ChangingAnyCounterChangesTheNumber)
{
	const real_t reference = standard_normal(7, 3, 11, 2);

	EXPECT_NE(standard_normal(8, 3, 11, 2), reference);
	EXPECT_NE(standard_normal(7, 4, 11, 2), reference);
	EXPECT_NE(standard_normal(7, 3, 12, 2), reference);
	EXPECT_NE(standard_normal(7, 3, 11, 1), reference);

	// swapping two counters is not the same as keeping them
	EXPECT_NE(standard_normal(7, 11, 3, 2), reference);
}

TEST(CounterRandomTest, NearbyCountersDoNotCollide)
{
	std::set<std::uint64_t> seen;
	for (std::uint64_t step = 0; step < 40; ++step)
	{
		for (std::uint64_t agent = 0; agent < 40; ++agent)
		{
			for (std::uint64_t component = 0; component < 3; ++component)
			{
				for (std::uint64_t stream = 0; stream < 2; ++stream)
				{
					seen.insert(random_bits(0, step, agent, component, stream));
				}
			}
		}
	}

	EXPECT_EQ(seen.size(), 40U * 40U * 3U * 2U);
}

TEST(CounterRandomTest, UnitIntervalNeverReturnsZero)
{
	EXPECT_GT(to_unit_interval(0), 0.0);
	EXPECT_DOUBLE_EQ(to_unit_interval(~std::uint64_t { 0 }), 1.0);
	EXPECT_TRUE(std::isfinite(std::log(to_unit_interval(0))));
}

TEST(CounterRandomTest, NumbersAreStandardNormalAcrossSubAgents)
{
	std::vector<real_t> samples(samples_count);
	for (index_t agent = 0; agent < samples_count; ++agent)
	{
		samples[agent] = standard_normal(42, 5, agent, 0);
	}

	const moments m = moments_of(samples);
	EXPECT_NEAR(m.mean, 0.0, five_sigma);
	EXPECT_NEAR(m.variance, 1.0, 0.02);
	EXPECT_NEAR(m.skewness, 0.0, 0.03);
	EXPECT_NEAR(m.kurtosis, 3.0, 0.06);
}

TEST(CounterRandomTest, NumbersAreStandardNormalAcrossSteps)
{
	// one sub-agent followed in time, the way a trajectory sees the numbers
	std::vector<real_t> samples(samples_count);
	for (index_t step = 0; step < samples_count; ++step)
	{
		samples[step] = standard_normal(42, step, 17, 1);
	}

	const moments m = moments_of(samples);
	EXPECT_NEAR(m.mean, 0.0, five_sigma);
	EXPECT_NEAR(m.variance, 1.0, 0.02);
	EXPECT_NEAR(m.skewness, 0.0, 0.03);
	EXPECT_NEAR(m.kurtosis, 3.0, 0.06);
}

TEST(CounterRandomTest, NumbersAreUncorrelated)
{
	std::vector<real_t> reference(samples_count);
	std::vector<real_t> next_step(samples_count);
	std::vector<real_t> next_agent(samples_count);
	std::vector<real_t> other_component(samples_count);
	std::vector<real_t> other_seed(samples_count);

	for (index_t i = 0; i < samples_count; ++i)
	{
		reference[i] = standard_normal(1, 100, i, 0);
		next_step[i] = standard_normal(1, 101, i, 0);
		next_agent[i] = standard_normal(1, 100, i + 1, 0);
		other_component[i] = standard_normal(1, 100, i, 1);
		other_seed[i] = standard_normal(2, 100, i, 0);
	}

	EXPECT_NEAR(correlation(reference, next_step), 0.0, five_sigma);
	EXPECT_NEAR(correlation(reference, next_agent), 0.0, five_sigma);
	EXPECT_NEAR(correlation(reference, other_component), 0.0, five_sigma);
	EXPECT_NEAR(correlation(reference, other_seed), 0.0, five_sigma);
}

} // namespace physicore::mechanics::micromechanics::kernels::openmp_solver::tests
