#include <gtest/gtest.h>

#include "micromechanics/cell_data.h"

namespace physicore::mechanics::micromechanics::tests {

TEST(CellDataTest, ConstructorSetsDimensions)
{
	const cell_data data(2);

	EXPECT_EQ(data.dims, index_t { 2 });
	EXPECT_EQ(data.cells_count, index_t { 0 });
	EXPECT_TRUE(data.positions.empty());
}

TEST(CellDataTest, AddInitializesAllFields)
{
	cell_data data(3);
	data.add();

	EXPECT_EQ(data.cells_count, index_t { 1 });
	EXPECT_EQ(data.is_movable, std::vector<std::uint8_t>({ 1 }));
	EXPECT_EQ(data.is_motile, std::vector<std::uint8_t>({ 0 }));
	EXPECT_EQ(data.migration_speeds, std::vector<real_t>({ 0.0 }));
	EXPECT_EQ(data.migration_biases, std::vector<real_t>({ 0.0 }));
	EXPECT_EQ(data.persistence_times, std::vector<real_t>({ 0.0 }));
	EXPECT_EQ(data.agent_counts, std::vector<index_t>({ 0 }));
	EXPECT_EQ(data.pressures, std::vector<real_t>({ 0.0 }));

	EXPECT_EQ(data.migration_bias_directions, std::vector<real_t>(3, 0.0));
	EXPECT_EQ(data.motility_directions, std::vector<real_t>(3, 0.0));
	EXPECT_EQ(data.positions, std::vector<real_t>(3, 0.0));
	EXPECT_EQ(data.velocities, std::vector<real_t>(3, 0.0));
	ASSERT_EQ(data.neighbor_cells.size(), 1U);
	EXPECT_TRUE(data.neighbor_cells[0].empty());
}

TEST(CellDataTest, RemoveMovesLastCellIntoRemovedSlot)
{
	cell_data data(2);
	data.add();
	data.add();
	data.add();

	data.is_movable = { 0, 1, 0 };
	data.is_motile = { 1, 0, 1 };
	data.migration_speeds = { 0.1, 0.2, 0.3 };
	data.migration_biases = { 0.4, 0.5, 0.6 };
	data.persistence_times = { 1.0, 2.0, 3.0 };
	data.agent_counts = { 10, 20, 30 };
	data.pressures = { 4.0, 5.0, 6.0 };
	data.migration_bias_directions = { 1.0, 1.1, 2.0, 2.1, 3.0, 3.1 };
	data.motility_directions = { 4.0, 4.1, 5.0, 5.1, 6.0, 6.1 };
	data.positions = { 7.0, 7.1, 8.0, 8.1, 9.0, 9.1 };
	data.velocities = { 10.0, 10.1, 11.0, 11.1, 12.0, 12.1 };
	data.neighbor_cells = { { 1 }, { 0, 2 }, { 42 } };

	data.remove_at(1);

	EXPECT_EQ(data.cells_count, index_t { 2 });
	EXPECT_EQ(data.is_movable, std::vector<std::uint8_t>({ 0, 0 }));
	EXPECT_EQ(data.is_motile, std::vector<std::uint8_t>({ 1, 1 }));
	EXPECT_EQ(data.migration_speeds, std::vector<real_t>({ 0.1, 0.3 }));
	EXPECT_EQ(data.migration_biases, std::vector<real_t>({ 0.4, 0.6 }));
	EXPECT_EQ(data.persistence_times, std::vector<real_t>({ 1.0, 3.0 }));
	EXPECT_EQ(data.agent_counts, std::vector<index_t>({ 10, 30 }));
	EXPECT_EQ(data.pressures, std::vector<real_t>({ 4.0, 6.0 }));
	EXPECT_EQ(data.migration_bias_directions, std::vector<real_t>({ 1.0, 1.1, 3.0, 3.1 }));
	EXPECT_EQ(data.motility_directions, std::vector<real_t>({ 4.0, 4.1, 6.0, 6.1 }));
	EXPECT_EQ(data.positions, std::vector<real_t>({ 7.0, 7.1, 9.0, 9.1 }));
	EXPECT_EQ(data.velocities, std::vector<real_t>({ 10.0, 10.1, 12.0, 12.1 }));
	EXPECT_EQ(data.neighbor_cells, std::vector<std::vector<index_t>>({ { 1 }, { 42 } }));
}

} // namespace physicore::mechanics::micromechanics::tests
