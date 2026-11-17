#include <array>
#include <cmath>

#include <gtest/gtest.h>

#include "micromechanics/pair_potential.h"
#include "micromechanics/simulation_parameters.h"

namespace physicore::mechanics::micromechanics::tests {

namespace {
// Deliberately not unit values, so a wrong power of a parameter cannot cancel out.
constexpr pair_potential potential {
	.well_depth = 3.0, .equilibrium_distance = 2.0, .steepness = 1.5, .cutoff_distance = 5.0
};

// Radial force (positive = repulsive) at distance r
real_t radial_force(const pair_potential& p, real_t r) { return force_coefficient(p, r * r) * r; }

// Distance at which t = exp(steepness * (1 - r^2 / r_eq^2)) takes the given value
real_t distance_for_t(const pair_potential& p, real_t t)
{
	return p.equilibrium_distance * std::sqrt(1.0 - std::log(t) / p.steepness);
}

real_t squared(real_t x) { return x * x; }
} // namespace

TEST(PairPotentialTest, MinimumIsAtEquilibriumDistanceWithWellDepth)
{
	const real_t r2 = squared(potential.equilibrium_distance);

	EXPECT_DOUBLE_EQ(potential_energy(potential, r2), -potential.well_depth);
	EXPECT_DOUBLE_EQ(force_coefficient(potential, r2), 0.0);

	// a minimum, not just a stationary point
	EXPECT_GT(potential_energy(potential, squared(potential.equilibrium_distance - 0.01)), -potential.well_depth);
	EXPECT_GT(potential_energy(potential, squared(potential.equilibrium_distance + 0.01)), -potential.well_depth);
}

TEST(PairPotentialTest, HandComputedValues)
{
	// 4 * well_depth * steepness / r_eq^2 = 4 * 3 * 1.5 / 4
	const real_t prefactor = 4.5;

	// r = 0: t = exp(steepness)
	const real_t t0 = std::exp(1.5);
	EXPECT_DOUBLE_EQ(potential_energy(potential, 0.0), 3.0 * t0 * (t0 - 2.0));
	EXPECT_DOUBLE_EQ(force_coefficient(potential, 0.0), prefactor * t0 * (t0 - 1.0));

	// t = 2: V = u0 * (4 - 4) = 0, c = prefactor * 2 * 1
	const real_t r_zero_energy = distance_for_t(potential, 2.0);
	EXPECT_NEAR(potential_energy(potential, squared(r_zero_energy)), 0.0, 1e-12);
	EXPECT_NEAR(force_coefficient(potential, squared(r_zero_energy)), 9.0, 1e-12);

	// t = 1/2: V = u0 * (1/4 - 1) = -2.25, c = prefactor * (1/2) * (-1/2) = -1.125
	const real_t r_half = distance_for_t(potential, 0.5);
	EXPECT_NEAR(potential_energy(potential, squared(r_half)), -2.25, 1e-12);
	EXPECT_NEAR(force_coefficient(potential, squared(r_half)), -1.125, 1e-12);
}

TEST(PairPotentialTest, RepulsiveBelowAndAttractiveAboveEquilibrium)
{
	for (const real_t r : { 0.0, 0.5, 1.0, 1.9, 1.999 })
	{
		EXPECT_GT(force_coefficient(potential, r * r), 0.0) << "r = " << r;
	}
	for (const real_t r : { 2.001, 2.5, 3.0, 4.0, 4.999 })
	{
		EXPECT_LT(force_coefficient(potential, r * r), 0.0) << "r = " << r;
	}
}

TEST(PairPotentialTest, ForceIsMinusDerivativeOfEnergy)
{
	const real_t h = 1e-5;
	for (const real_t r : { 0.3, 1.0, 1.7, 2.0, 2.4, 3.1, 4.5 })
	{
		const real_t numerical =
			-(potential_energy(potential, squared(r + h)) - potential_energy(potential, squared(r - h))) / (2.0 * h);

		EXPECT_NEAR(radial_force(potential, r), numerical, 1e-7 * (1.0 + std::abs(numerical))) << "r = " << r;
	}
}

TEST(PairPotentialTest, StiffnessAtEquilibrium)
{
	// V''(r_eq) = 8 * well_depth * steepness^2 / r_eq^2 = 8 * 3 * 2.25 / 4
	const real_t expected_stiffness = 13.5;

	const real_t h = 1e-5;
	const real_t r_eq = potential.equilibrium_distance;
	const real_t numerical = -(radial_force(potential, r_eq + h) - radial_force(potential, r_eq - h)) / (2.0 * h);

	EXPECT_NEAR(numerical, expected_stiffness, 1e-6);
}

TEST(PairPotentialTest, NoInteractionFromCutoffOn)
{
	const real_t cutoff2 = squared(potential.cutoff_distance);

	EXPECT_NE(force_coefficient(potential, squared(potential.cutoff_distance - 1e-6)), 0.0);
	EXPECT_NE(potential_energy(potential, squared(potential.cutoff_distance - 1e-6)), 0.0);

	EXPECT_DOUBLE_EQ(force_coefficient(potential, cutoff2), 0.0);
	EXPECT_DOUBLE_EQ(potential_energy(potential, cutoff2), 0.0);
	EXPECT_DOUBLE_EQ(force_coefficient(potential, 100.0 * cutoff2), 0.0);
	EXPECT_DOUBLE_EQ(potential_energy(potential, 100.0 * cutoff2), 0.0);
}

TEST(PairPotentialTest, FiniteAtZeroDistance)
{
	EXPECT_TRUE(std::isfinite(potential_energy(potential, 0.0)));
	EXPECT_TRUE(std::isfinite(force_coefficient(potential, 0.0)));
	// two coincident elements feel no force: the direction is undefined and c * (x_i - x_j) is zero
	EXPECT_DOUBLE_EQ(radial_force(potential, 0.0), 0.0);
}

TEST(PairPotentialTest, WellDepthScalesEnergyAndForceLinearly)
{
	pair_potential deeper = potential;
	deeper.well_depth = 2.0 * potential.well_depth;

	for (const real_t r : { 0.5, 2.0, 3.0 })
	{
		EXPECT_DOUBLE_EQ(potential_energy(deeper, r * r), 2.0 * potential_energy(potential, r * r));
		EXPECT_DOUBLE_EQ(force_coefficient(deeper, r * r), 2.0 * force_coefficient(potential, r * r));
	}
}

TEST(PairPotentialTest, SignConventionAndNewtonsThirdLaw)
{
	// elements closer than the equilibrium distance: i sits at larger x than j
	const std::array<real_t, 3> xi = { 1.0, 0.5, 0.0 };
	const std::array<real_t, 3> xj = { 0.0, 0.0, 0.0 };

	const real_t c = force_coefficient(potential, 1.25);
	ASSERT_GT(c, 0.0);

	for (std::size_t d = 0; d < 3; ++d)
	{
		const real_t force_on_i = c * (xi.at(d) - xj.at(d));
		const real_t force_on_j = c * (xj.at(d) - xi.at(d));

		EXPECT_DOUBLE_EQ(force_on_i + force_on_j, 0.0);
		// repulsion moves i further away from j along every axis where they differ
		EXPECT_GE(force_on_i * (xi.at(d) - xj.at(d)), 0.0);
	}
}

TEST(PairPotentialTest, SameSeparationGivesSameForceIn2DAnd3D)
{
	// the same distance 2.5, once in the plane and once along the space diagonal of a box
	const real_t in_plane = squared(1.5) + squared(2.0);
	const real_t in_space = squared(1.5) + squared(1.2) + squared(1.6);

	EXPECT_DOUBLE_EQ(in_plane, in_space);
	EXPECT_DOUBLE_EQ(force_coefficient(potential, in_plane), force_coefficient(potential, in_space));
}

TEST(SimulationParametersTest, IntraAndInterPotentialsShareShapeAndDifferInDepth)
{
	simulation_parameters params;
	params.equilibrium_distance = 2.0;
	params.steepness = 1.5;
	params.intra_well_depth = 3.0;
	params.inter_well_depth = 0.75;
	params.cutoff_factor = 2.5;

	const pair_potential intra = params.intra_potential();
	const pair_potential inter = params.inter_potential();

	EXPECT_DOUBLE_EQ(intra.well_depth, 3.0);
	EXPECT_DOUBLE_EQ(inter.well_depth, 0.75);
	for (const pair_potential& p : { intra, inter })
	{
		EXPECT_DOUBLE_EQ(p.equilibrium_distance, 2.0);
		EXPECT_DOUBLE_EQ(p.steepness, 1.5);
		EXPECT_DOUBLE_EQ(p.cutoff_distance, 5.0);
	}
}

TEST(SimulationParametersTest, DefaultsKeepCellsCohesive)
{
	const simulation_parameters params;

	EXPECT_GT(params.intra_well_depth, params.inter_well_depth);
	EXPECT_DOUBLE_EQ(params.steepness, 2.0);
	EXPECT_DOUBLE_EQ(params.cutoff_factor, 2.5);
}

} // namespace physicore::mechanics::micromechanics::tests
