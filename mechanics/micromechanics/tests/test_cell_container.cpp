#include <memory>

#include <gtest/gtest.h>

#include "micromechanics/cell_container.h"

namespace physicore::mechanics::micromechanics::tests {

namespace {
cell_container make_cell_container() { return cell_container(std::make_unique<cell_data>(3)); }
} // namespace

TEST(CellContainerTest, CreateAndAccessCells)
{
	cell_container container = make_cell_container();

	auto* cell0 = container.create();
	auto* cell1 = container.create();

	ASSERT_NE(cell0, nullptr);
	ASSERT_NE(cell1, nullptr);
	EXPECT_EQ(container.size(), 2U);

	cell0->velocity()[0] = 1.5;
	cell1->velocity()[0] = 2.5;

	EXPECT_EQ(container.get_agent_at(0), cell0);
	EXPECT_EQ(container.get_agent_at(1), cell1);
	EXPECT_DOUBLE_EQ(container.get_agent_at(0)->velocity()[0], 1.5);
	EXPECT_DOUBLE_EQ(container.get_agent_at(1)->velocity()[0], 2.5);
}

TEST(CellContainerTest, RemoveMovesLastCellAndItsData)
{
	cell_container container = make_cell_container();

	auto* cell0 = container.create();
	auto* removed_cell = container.create();
	auto* cell2 = container.create();

	cell0->velocity()[0] = 1.0;
	removed_cell->velocity()[0] = 2.0;
	cell2->velocity()[0] = 3.0;
	cell2->position()[0] = 4.0;

	container.remove_agent(removed_cell);

	ASSERT_EQ(container.size(), 2U);
	EXPECT_EQ(container.get_agent_at(0), cell0);
	EXPECT_EQ(container.get_agent_at(1), cell2);
	EXPECT_DOUBLE_EQ(container.get_agent_at(1)->velocity()[0], 3.0);
	EXPECT_DOUBLE_EQ(container.get_agent_at(1)->position()[0], 4.0);
}

} // namespace physicore::mechanics::micromechanics::tests
