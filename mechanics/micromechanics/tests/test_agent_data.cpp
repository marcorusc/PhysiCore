#include <vector>

#include <common/base_agent_data.h>
#include <gtest/gtest.h>

#include "micromechanics/agent_data.h"
#include "micromechanics/agent_interface.h"

namespace physicore::mechanics::micromechanics::tests {

namespace {
void expect_sizes(const agent_data& data, index_t count, index_t dims)
{
	EXPECT_EQ(data.agents_count, count);
	EXPECT_EQ(data.cell_indices.size(), count);
	EXPECT_EQ(data.velocities.size(), count * dims);
}

void sync_add(base_agent_data& base, agent_data& data)
{
	base.add();
	data.add();
}

void sync_remove_at(base_agent_data& base, agent_data& data, index_t position)
{
	base.remove_at(position);
	data.remove_at(position);
}

class AgentDataTest : public ::testing::TestWithParam<index_t>
{};
} // namespace

INSTANTIATE_TEST_SUITE_P(Dimensions, AgentDataTest, ::testing::Values(1, 2, 3));

TEST_P(AgentDataTest, AddGrowsAllArraysWithDefaults)
{
	const index_t dims = GetParam();
	base_agent_data base(dims);
	agent_data data(base);

	expect_sizes(data, 0, dims);

	sync_add(base, data);
	sync_add(base, data);

	expect_sizes(data, 2, dims);
	EXPECT_EQ(data.cell_indices, std::vector<index_t>({ no_cell, no_cell }));
	EXPECT_EQ(data.velocities, std::vector<real_t>(2 * dims, 0.0));
}

TEST_P(AgentDataTest, RemoveMovesLastAgentIntoRemovedSlot)
{
	const index_t dims = GetParam();
	base_agent_data base(dims);
	agent_data data(base);

	for (index_t i = 0; i < 3; ++i)
	{
		sync_add(base, data);
		data.cell_indices[i] = 10 + i;
		for (index_t d = 0; d < dims; ++d)
		{
			data.velocities[i * dims + d] = static_cast<real_t>(100 * i + d);
		}
	}

	sync_remove_at(base, data, 0);

	expect_sizes(data, 2, dims);
	EXPECT_EQ(data.cell_indices, std::vector<index_t>({ 12, 11 }));
	for (index_t d = 0; d < dims; ++d)
	{
		EXPECT_DOUBLE_EQ(data.velocities[d], static_cast<real_t>(200 + d));
		EXPECT_DOUBLE_EQ(data.velocities[dims + d], static_cast<real_t>(100 + d));
	}

	// removing the last agent moves nothing
	sync_remove_at(base, data, 1);

	expect_sizes(data, 1, dims);
	EXPECT_EQ(data.cell_indices, std::vector<index_t>({ 12 }));

	sync_remove_at(base, data, 0);

	expect_sizes(data, 0, dims);
}

} // namespace physicore::mechanics::micromechanics::tests
