#include <cmath>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "micromechanics/cell_type.h"
#include "micromechanics/environment.h"

// A cell type says what the elements of a cell are like and how its cells stick to others. The
// interaction_table works out, once, the interaction of every pair of types. These tests check the
// rule for each of the three kinds of pairs, the weighting of the attraction by affinities, and that
// the environment creates cells of a type.

namespace physicore::mechanics::micromechanics::tests {

namespace {
// Two types with different element sizes, well depths and affinities
std::vector<cell_type> two_types()
{
	cell_type small;
	small.name = "small";
	small.equilibrium_distance = 1.0;
	small.intra_well_depth = 2.0;
	small.inter_well_depth = 0.5;
	small.adhesion_affinities = { 0.9, 0.4 };

	cell_type large;
	large.name = "large";
	large.equilibrium_distance = 3.0;
	large.intra_well_depth = 5.0;
	large.inter_well_depth = 2.0;
	large.adhesion_affinities = { 0.1, 0.6 };

	return { small, large };
}

void expect_same(const pair_interaction& a, const pair_interaction& b)
{
	EXPECT_DOUBLE_EQ(a.potential.well_depth, b.potential.well_depth);
	EXPECT_DOUBLE_EQ(a.potential.equilibrium_distance, b.potential.equilibrium_distance);
	EXPECT_DOUBLE_EQ(a.potential.steepness, b.potential.steepness);
	EXPECT_DOUBLE_EQ(a.potential.cutoff_distance, b.potential.cutoff_distance);
	EXPECT_DOUBLE_EQ(a.attraction_scale, b.attraction_scale);
}

environment make_environment(index_t dims = 3)
{
	return { { dims, { -100, -100, -100 }, { 100, 100, 100 }, { 20, 20, 20 } }, 0.01 };
}
} // namespace

// ------------------------------------------------------------------------------------ the type

TEST(CellTypeTest, DefaultsAreUnitValuesWithTheInterCellularWellHalfAsDeep)
{
	const cell_type type;

	EXPECT_DOUBLE_EQ(type.drag_coefficient, 1.0);
	EXPECT_DOUBLE_EQ(type.equilibrium_distance, 1.0);
	EXPECT_DOUBLE_EQ(type.intra_well_depth, 1.0);
	EXPECT_DOUBLE_EQ(type.inter_well_depth, 0.5);
	EXPECT_TRUE(type.adhesion_affinities.empty());
	EXPECT_FALSE(type.is_motile);
}

TEST(CellTypeTest, AffinityWithoutAnEntryIsOne)
{
	cell_type type;
	type.adhesion_affinities = { 0.3 };

	EXPECT_DOUBLE_EQ(type.affinity_to(0), 0.3);
	EXPECT_DOUBLE_EQ(type.affinity_to(1), 1.0);
	EXPECT_DOUBLE_EQ(type.affinity_to(7), 1.0);
}

// ------------------------------------------------------------------------ the weighted potential

TEST(PairInteractionTest, RepulsionIsUntouchedAndAttractionIsScaled)
{
	const pair_potential potential {
		.well_depth = 2.0, .equilibrium_distance = 1.5, .steepness = 2.0, .cutoff_distance = 3.75
	};
	const pair_interaction weak { .potential = potential, .attraction_scale = 0.25 };

	for (const real_t distance : { 0.3, 0.9, 1.4 })
	{
		EXPECT_DOUBLE_EQ(force_coefficient(weak, distance * distance),
						 force_coefficient(potential, distance * distance))
			<< "distance " << distance;
		EXPECT_GT(force_coefficient(weak, distance * distance), 0.0);
	}
	for (const real_t distance : { 1.6, 2.0, 3.0 })
	{
		EXPECT_DOUBLE_EQ(force_coefficient(weak, distance * distance),
						 0.25 * force_coefficient(potential, distance * distance))
			<< "distance " << distance;
		EXPECT_LT(force_coefficient(weak, distance * distance), 0.0);
	}

	// no force at the equilibrium distance, whatever the scale, and none beyond the cutoff
	EXPECT_NEAR(force_coefficient(weak, 1.5 * 1.5), 0.0, 1e-15);
	EXPECT_DOUBLE_EQ(force_coefficient(weak, 4.0 * 4.0), 0.0);
}

TEST(PairInteractionTest, ScaleOneIsThePlainPotentialAndScaleZeroOnlyRepels)
{
	const pair_potential potential {
		.well_depth = 1.0, .equilibrium_distance = 1.0, .steepness = 2.0, .cutoff_distance = 2.5
	};
	const pair_interaction full { .potential = potential, .attraction_scale = 1.0 };
	const pair_interaction none { .potential = potential, .attraction_scale = 0.0 };

	for (int k = 1; k < 60; ++k)
	{
		const real_t squared = std::pow(0.05 * k, 2);
		EXPECT_DOUBLE_EQ(force_coefficient(full, squared), force_coefficient(potential, squared));
		EXPECT_DOUBLE_EQ(potential_energy(full, squared), potential_energy(potential, squared));
		EXPECT_GE(force_coefficient(none, squared), 0.0) << "distance " << 0.05 * k;
		EXPECT_GE(potential_energy(none, squared), 0.0) << "distance " << 0.05 * k;
	}
}

TEST(PairInteractionTest, ForceIsMinusTheDerivativeOfTheEnergyOnBothSidesOfTheMinimum)
{
	const pair_potential potential {
		.well_depth = 2.0, .equilibrium_distance = 1.5, .steepness = 2.0, .cutoff_distance = 3.75
	};
	const pair_interaction weak { .potential = potential, .attraction_scale = 0.3 };
	const real_t step = 1e-6;

	for (const real_t distance : { 0.8, 1.2, 1.45, 1.55, 2.0, 2.8 })
	{
		const real_t derivative = (potential_energy(weak, std::pow(distance + step, 2))
								   - potential_energy(weak, std::pow(distance - step, 2)))
								  / (2.0 * step);
		// the radial force is coefficient * distance
		EXPECT_NEAR(force_coefficient(weak, distance * distance) * distance, -derivative, 1e-6)
			<< "distance " << distance;
	}

	// the energy has no jump at the equilibrium distance: the bottom of the well is scale * (-depth)
	EXPECT_NEAR(potential_energy(weak, std::pow(1.5 - 1e-9, 2)), -0.3 * 2.0, 1e-7);
	EXPECT_NEAR(potential_energy(weak, std::pow(1.5 + 1e-9, 2)), -0.3 * 2.0, 1e-7);
}

// --------------------------------------------------------------------------------------- the table

TEST(InteractionTableTest, ElementsOfOneCellUseTheIntraCellularWellOfItsType)
{
	const auto types = two_types();
	const simulation_parameters params;
	const interaction_table table(types, params);

	ASSERT_EQ(table.types_count(), index_t { 2 });
	expect_same(table.intra(0), { .potential = intra_potential(types[0], params), .attraction_scale = 1.0 });
	expect_same(table.intra(1), { .potential = intra_potential(types[1], params), .attraction_scale = 1.0 });
	EXPECT_DOUBLE_EQ(table.intra(1).potential.well_depth, 5.0);
	EXPECT_DOUBLE_EQ(table.intra(1).potential.cutoff_distance, 2.5 * 3.0);
}

TEST(InteractionTableTest, CellsOfTheSameTypeUseItsInterCellularWellWeightedByItsAffinityToItself)
{
	const auto types = two_types();
	const simulation_parameters params;
	const interaction_table table(types, params);

	expect_same(table.inter(0, 0), { .potential = inter_potential(types[0], params), .attraction_scale = 0.9 });
	expect_same(table.inter(1, 1), { .potential = inter_potential(types[1], params), .attraction_scale = 0.6 });
}

TEST(InteractionTableTest, CellsOfDifferentTypesUseAWellMixedFromBoth)
{
	const auto types = two_types();
	const simulation_parameters params;
	const interaction_table table(types, params);

	const pair_interaction& mixed = table.inter(0, 1);
	// the sum of the two element radii: (1 + 3) / 2
	EXPECT_DOUBLE_EQ(mixed.potential.equilibrium_distance, 2.0);
	// geometric mean of the inter-cellular well depths: sqrt(0.5 * 2)
	EXPECT_DOUBLE_EQ(mixed.potential.well_depth, 1.0);
	// geometric mean of the two affinities: sqrt(0.4 * 0.1)
	EXPECT_DOUBLE_EQ(mixed.attraction_scale, 0.2);
	EXPECT_DOUBLE_EQ(mixed.potential.cutoff_distance, 2.5 * 2.0);
	EXPECT_DOUBLE_EQ(mixed.potential.steepness, 2.0);
}

TEST(InteractionTableTest, InteractionIsTheSameSeenFromEitherTypeEvenWithUnequalAffinities)
{
	// the affinity of small to large is 0.4 and of large to small 0.1: the table has one value for both
	const auto types = two_types();
	const interaction_table table(types, simulation_parameters {});

	expect_same(table.inter(0, 1), table.inter(1, 0));
}

TEST(InteractionTableTest, LargestCutoffCoversEveryPair)
{
	const auto types = two_types();
	const interaction_table table(types, simulation_parameters {});

	// the widest pair is two elements of the large type
	EXPECT_DOUBLE_EQ(table.largest_cutoff(), 2.5 * 3.0);
}

TEST(InteractionTableTest, RejectsTypesThatMakeNoSense)
{
	const simulation_parameters params;
	EXPECT_THROW(interaction_table(std::vector<cell_type> {}, params), std::invalid_argument);

	const auto broken = [&](auto&& change) {
		std::vector<cell_type> types = two_types();
		change(types[1]);
		return [types, params] { const interaction_table table(types, params); };
	};

	EXPECT_THROW(broken([](cell_type& t) { t.drag_coefficient = 0.0; })(), std::invalid_argument);
	EXPECT_THROW(broken([](cell_type& t) { t.equilibrium_distance = -1.0; })(), std::invalid_argument);
	EXPECT_THROW(broken([](cell_type& t) { t.intra_well_depth = 0.0; })(), std::invalid_argument);
	EXPECT_THROW(broken([](cell_type& t) { t.inter_well_depth = -0.1; })(), std::invalid_argument);
	EXPECT_THROW(broken([](cell_type& t) { t.adhesion_affinities = { 1.0, 1.5 }; })(), std::invalid_argument);
	EXPECT_THROW(broken([](cell_type& t) { t.adhesion_affinities = { -0.2 }; })(), std::invalid_argument);

	// cells that do not attract or even feel each other at all are allowed
	EXPECT_NO_THROW(broken([](cell_type& t) { t.inter_well_depth = 0.0; })());
	EXPECT_NO_THROW(broken([](cell_type& t) { t.adhesion_affinities = { 0.0, 0.0 }; })());

	simulation_parameters flat;
	flat.steepness = 0.0;
	EXPECT_THROW(interaction_table(two_types(), flat), std::invalid_argument);
}

// --------------------------------------------------------------------------------- the environment

TEST(CellTypeEnvironmentTest, EnvironmentStartsWithOneDefaultTypeAndCanGetMore)
{
	environment env = make_environment();

	ASSERT_EQ(env.cell_types.size(), 1U);
	EXPECT_DOUBLE_EQ(env.cell_types[0].equilibrium_distance, 1.0);

	EXPECT_EQ(env.add_cell_type(two_types()[1]), index_t { 1 });
	EXPECT_EQ(env.add_cell_type(two_types()[0]), index_t { 2 });
	EXPECT_EQ(env.cell_types[1].name, "large");
}

TEST(CellTypeEnvironmentTest, NewCellTakesTheTypeItIsCreatedWith)
{
	environment env = make_environment();
	const index_t large_index = env.add_cell_type(two_types()[1]);

	EXPECT_EQ(env.create_cell()->type_index(), index_t { 0 });
	EXPECT_EQ(env.create_cell(large_index)->type_index(), large_index);

	EXPECT_THROW(env.create_cell(5), std::out_of_range);
}

TEST(CellTypeEnvironmentTest, NewCellTakesTheMotilityOfItsType)
{
	environment env = make_environment();
	cell_type crawler;
	crawler.is_motile = true;
	crawler.migration_speed = 0.7;
	crawler.persistence_time = 3.0;
	crawler.migration_bias = 0.25;
	const index_t crawler_index = env.add_cell_type(crawler);

	auto* plain = env.create_cell();
	auto* motile = env.create_cell(crawler_index);

	EXPECT_EQ(plain->is_motile(), 0);

	EXPECT_EQ(motile->is_motile(), 1);
	EXPECT_DOUBLE_EQ(motile->migration_speed(), 0.7);
	EXPECT_DOUBLE_EQ(motile->persistence_time(), 3.0);
	EXPECT_DOUBLE_EQ(motile->migration_bias(), 0.25);
}

} // namespace physicore::mechanics::micromechanics::tests
