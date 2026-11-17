#include <memory>

#include <common/base_agent_container.h>
#include <gtest/gtest.h>

#include "micromechanics/agent_container.h"

namespace physicore::mechanics::micromechanics::tests {

namespace {
std::shared_ptr<agent_container> make_agent_container()
{
	auto base_data = std::make_unique<base_agent_data>(3);
	auto data = std::make_unique<agent_data>(*base_data);
	return std::make_shared<agent_container>(std::move(base_data), std::move(data));
}
} // namespace

TEST(AgentContainerTest, CreateAndAccessAgents)
{
	auto container = make_agent_container();

	auto* agent0 = container->create();
	auto* agent1 = container->create();

	ASSERT_NE(agent0, nullptr);
	ASSERT_NE(agent1, nullptr);
	EXPECT_EQ(container->size(), 2U);

	agent0->cell_index() = 0;
	agent1->cell_index() = 1;
	agent1->velocity()[2] = 5.0;

	EXPECT_EQ(container->get_agent_at(0), agent0);
	EXPECT_EQ(container->get_agent_at(1), agent1);
	EXPECT_EQ(container->get_agent_at(0)->cell_index(), index_t { 0 });
	EXPECT_EQ(container->get_agent_at(1)->cell_index(), index_t { 1 });
	EXPECT_DOUBLE_EQ(container->get_agent_at(1)->velocity()[2], 5.0);
}

TEST(AgentContainerTest, RemoveMovesLastAgentAndItsData)
{
	auto container = make_agent_container();

	auto* agent0 = container->create();
	auto* removed_agent = container->create();
	auto* agent2 = container->create();

	agent0->cell_index() = 0;
	removed_agent->cell_index() = 1;
	agent2->cell_index() = 2;
	agent2->position()[0] = 4.0;
	agent2->velocity()[1] = 6.0;

	container->remove_agent(removed_agent);

	ASSERT_EQ(container->size(), 2U);
	EXPECT_EQ(container->get_agent_at(0), agent0);
	EXPECT_EQ(container->get_agent_at(1), agent2);
	EXPECT_EQ(agent2->cell_index(), index_t { 2 });
	EXPECT_DOUBLE_EQ(agent2->position()[0], 4.0);
	EXPECT_DOUBLE_EQ(agent2->velocity()[1], 6.0);
}

TEST(AgentContainerTest, ViewsShareOneSetOfAgents)
{
	auto container = make_agent_container();

	const agent_container_ptr mechanics_view = container;
	physicore::generic_agent_interface_container<base_agent_interface>& base_view = *container;

	auto* created = mechanics_view->create();
	created->position()[1] = 3.0;

	ASSERT_EQ(base_view.size(), 1U);
	EXPECT_DOUBLE_EQ(base_view.get_agent_at(0)->position()[1], 3.0);

	base_view.create();
	EXPECT_EQ(mechanics_view->size(), 2U);
	EXPECT_EQ(mechanics_view->get_agent_at(1)->cell_index(), no_cell);

	base_view.remove_at(0);
	EXPECT_EQ(mechanics_view->size(), 1U);
}

} // namespace physicore::mechanics::micromechanics::tests
