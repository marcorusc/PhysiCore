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
	data.agent_counts[0] = 9;
	data.neighbor_cells[0] = { 11, 12 };

	EXPECT_DOUBLE_EQ(data.positions[0], 1.0);
	EXPECT_DOUBLE_EQ(data.velocities[1], 2.0);
	EXPECT_EQ(data.is_movable[0], std::uint8_t { 0 });
	EXPECT_EQ(data.agent_counts[0], index_t { 9 });
	ASSERT_EQ(test_cell.neighbors().size(), 2U);
	EXPECT_EQ(test_cell.neighbors()[0], index_t { 11 });
	EXPECT_EQ(test_cell.neighbors()[1], index_t { 12 });
	EXPECT_EQ(test_cell.agent_count(), index_t { 9 });
}

} // namespace physicore::mechanics::micromechanics::tests
