#include <cmath>
#include <numbers>
#include <stdexcept>

#include <gtest/gtest.h>

#include "micromechanics/cell_parameters.h"

namespace physicore::mechanics::micromechanics::tests {

namespace {
// The cell of Sandersius and Newman (2008), section 4.1, in micrometres, seconds and piconewtons
cell_parameters paper_cell() { return { .radius = 10.0, .stiffness = 5000.0, .damping = 5000.0 }; }

real_t pair_stiffness(const simulation_parameters& params)
{
	return 8.0 * params.intra_well_depth * params.steepness * params.steepness
		   / (params.equilibrium_distance * params.equilibrium_distance);
}
} // namespace

TEST(CellParametersTest, DefaultsAreThoseOfThePaper)
{
	const cell_parameters cell;

	EXPECT_DOUBLE_EQ(cell.packing_density, std::numbers::pi / (3.0 * std::numbers::sqrt2));
	EXPECT_NEAR(cell.packing_density, 0.7405, 5e-5);
	EXPECT_DOUBLE_EQ(cell.scaling_correction, 0.75);
}

TEST(CellParametersTest, ReproducesTheEquilibriumDistanceQuotedInThePaper)
{
	// "r_eq = 2.85 um as computed from equation (8) with N = 256, R_cell = 10 um and p_3 = pi / (3 sqrt 2)"
	const simulation_parameters params = scale_to_elements(paper_cell(), 256);

	EXPECT_NEAR(params.equilibrium_distance, 2.85, 0.005);
}

TEST(CellParametersTest, HandComputedScalingForAThousandElements)
{
	// N = 1000, so N^(-1/3) = 0.1
	const simulation_parameters params = scale_to_elements(paper_cell(), 1000);

	// Eq. 8: 2 * 10 * (p / 1000)^(1/3) = 2 * cbrt(p)
	const real_t expected_distance = 2.0 * std::cbrt(std::numbers::pi / (3.0 * std::numbers::sqrt2));
	EXPECT_NEAR(params.equilibrium_distance, expected_distance, 1e-12);

	// Eq. 12: 5000 / 1000
	EXPECT_NEAR(params.drag_coefficient, 5.0, 1e-12);

	// Eq. 11: 5000 * 0.1 * (1 - 0.75 * 0.1) = 462.5, and it is the stiffness of the resulting potential
	EXPECT_NEAR(pair_stiffness(params), 462.5, 1e-9);

	// Eq. 5 with the default steepness of 2: u0 = 462.5 * r_eq^2 / 32
	EXPECT_NEAR(params.intra_well_depth, 462.5 * expected_distance * expected_distance / 32.0, 1e-9);
}

TEST(CellParametersTest, ElementsFillThePackedFractionOfTheCellVolume)
{
	// Eq. 7: N spheres of diameter r_eq occupy the fraction p_3 of the volume of the cell
	const cell_parameters cell = paper_cell();
	for (const index_t elements : { 100, 256, 1000, 4096 })
	{
		const simulation_parameters params = scale_to_elements(cell, elements);

		const real_t elements_volume = static_cast<real_t>(elements) * std::pow(0.5 * params.equilibrium_distance, 3);
		const real_t cell_volume = std::pow(cell.radius, 3);

		EXPECT_NEAR(elements_volume / cell_volume, cell.packing_density, 1e-12) << elements << " elements";
	}
}

TEST(CellParametersTest, MoreElementsAreSmallerSofterAndLessDamped)
{
	const simulation_parameters coarse = scale_to_elements(paper_cell(), 256);
	const simulation_parameters fine = scale_to_elements(paper_cell(), 2048);

	// eight times as many elements: half the spacing, an eighth of the drag
	EXPECT_NEAR(fine.equilibrium_distance, 0.5 * coarse.equilibrium_distance, 1e-12);
	EXPECT_NEAR(fine.drag_coefficient, coarse.drag_coefficient / 8.0, 1e-12);
	EXPECT_LT(pair_stiffness(fine), pair_stiffness(coarse));
	EXPECT_LT(fine.intra_well_depth, coarse.intra_well_depth);
}

TEST(CellParametersTest, ShapeOfThePotentialIsTakenOver)
{
	simulation_parameters shape;
	shape.steepness = 4.0;
	shape.cutoff_factor = 2.0;
	shape.intra_well_depth = 2.0;
	shape.inter_well_depth = 0.5;

	const simulation_parameters params = scale_to_elements(paper_cell(), 1000, shape);

	EXPECT_DOUBLE_EQ(params.steepness, 4.0);
	EXPECT_DOUBLE_EQ(params.cutoff_factor, 2.0);
	EXPECT_NEAR(params.inter_well_depth / params.intra_well_depth, 0.25, 1e-12);
	// the pair stiffness does not depend on the steepness, the well depth compensates
	EXPECT_NEAR(pair_stiffness(params), 462.5, 1e-9);
}

TEST(CellParametersTest, RejectsValuesThatMakeNoSense)
{
	EXPECT_THROW(scale_to_elements(paper_cell(), 0), std::invalid_argument);

	cell_parameters no_radius = paper_cell();
	no_radius.radius = 0.0;
	EXPECT_THROW(scale_to_elements(no_radius, 256), std::invalid_argument);

	cell_parameters negative_stiffness = paper_cell();
	negative_stiffness.stiffness = -1.0;
	EXPECT_THROW(scale_to_elements(negative_stiffness, 256), std::invalid_argument);

	// a default-constructed cell has no size or stiffness yet
	EXPECT_THROW(scale_to_elements(cell_parameters {}, 256), std::invalid_argument);

	// a single element with the default correction would get the stiffness 1 - 0.75 > 0, but a large
	// correction turns the stiffness of a coarse cell negative
	cell_parameters large_correction = paper_cell();
	large_correction.scaling_correction = 3.0;
	EXPECT_THROW(scale_to_elements(large_correction, 8), std::invalid_argument);
	EXPECT_NO_THROW(scale_to_elements(large_correction, 1000));
}

} // namespace physicore::mechanics::micromechanics::tests
