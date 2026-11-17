#include <common/base_agent_data.h>
#include <gtest/gtest.h>

#include "micromechanics/agent.h"

namespace physicore::mechanics::micromechanics::tests {

TEST(AgentTest, AccessorsReferenceAgentData)
{
	const index_t dims = 3;
	base_agent_data base(dims);
	agent_data data(base);

	for (index_t i = 0; i < 2; ++i)
	{
		base.add();
		data.add();
	}

	agent second(1, data);

	EXPECT_EQ(second.cell_index(), no_cell);
	ASSERT_EQ(second.position().size(), dims);
	ASSERT_EQ(second.velocity().size(), dims);
	ASSERT_EQ(second.previous_velocity().size(), dims);

	second.cell_index() = 7;
	for (index_t d = 0; d < dims; ++d)
	{
		second.position()[d] = static_cast<real_t>(1 + d);
		second.velocity()[d] = static_cast<real_t>(10 + d);
		second.previous_velocity()[d] = static_cast<real_t>(20 + d);
	}

	EXPECT_EQ(data.cell_indices[1], index_t { 7 });
	for (index_t d = 0; d < dims; ++d)
	{
		EXPECT_DOUBLE_EQ(base.positions[dims + d], static_cast<real_t>(1 + d));
		EXPECT_DOUBLE_EQ(data.velocities[dims + d], static_cast<real_t>(10 + d));
		EXPECT_DOUBLE_EQ(data.previous_velocities[dims + d], static_cast<real_t>(20 + d));
	}

	// the first agent is untouched
	EXPECT_EQ(data.cell_indices[0], no_cell);
	for (index_t d = 0; d < dims; ++d)
	{
		EXPECT_DOUBLE_EQ(base.positions[d], 0.0);
		EXPECT_DOUBLE_EQ(data.velocities[d], 0.0);
		EXPECT_DOUBLE_EQ(data.previous_velocities[d], 0.0);
	}
}

} // namespace physicore::mechanics::micromechanics::tests
