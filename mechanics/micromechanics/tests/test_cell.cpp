#include <gtest/gtest.h>

#include "micromechanics/cell.h"

namespace physicore::mechanics::micromechanics::tests {

class CellTest : public ::testing::Test
{
protected:
	cell_data data { 3 };

	void SetUp() override { data.add(); }
};

TEST_F(CellTest, AccessorsReferenceCellData)
{
	cell test_cell(0, data);

	test_cell.position()[0] = 1.0;
	test_cell.velocity()[1] = 2.0;
	test_cell.is_movable() = 0;
	test_cell.is_motile() = 1;
	test_cell.migration_speed() = 4.0;
	test_cell.migration_bias() = 0.5;
	test_cell.migration_bias_direction()[2] = 6.0;
	test_cell.persistence_time() = 7.0;
	test_cell.motility_direction()[0] = 8.0;
	data.agent_counts[0] = 9;
	data.pressures[0] = 10.0;
	data.neighbor_cells[0] = { 11, 12 };

	EXPECT_DOUBLE_EQ(data.positions[0], 1.0);
	EXPECT_DOUBLE_EQ(data.velocities[1], 2.0);
	EXPECT_EQ(data.is_movable[0], std::uint8_t { 0 });
	EXPECT_EQ(data.is_motile[0], std::uint8_t { 1 });
	EXPECT_DOUBLE_EQ(data.migration_speeds[0], 4.0);
	EXPECT_DOUBLE_EQ(data.migration_biases[0], 0.5);
	EXPECT_DOUBLE_EQ(data.migration_bias_directions[2], 6.0);
	EXPECT_DOUBLE_EQ(data.persistence_times[0], 7.0);
	EXPECT_DOUBLE_EQ(data.motility_directions[0], 8.0);
	EXPECT_EQ(data.agent_counts[0], index_t { 9 });
	EXPECT_DOUBLE_EQ(data.pressures[0], 10.0);
	ASSERT_EQ(test_cell.neighbors().size(), 2U);
	EXPECT_EQ(test_cell.neighbors()[0], index_t { 11 });
	EXPECT_EQ(test_cell.neighbors()[1], index_t { 12 });
	EXPECT_EQ(test_cell.agent_count(), index_t { 9 });
	EXPECT_DOUBLE_EQ(test_cell.pressure(), 10.0);
}

} // namespace physicore::mechanics::micromechanics::tests
