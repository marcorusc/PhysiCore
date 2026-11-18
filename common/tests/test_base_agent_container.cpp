#include <memory>
#include <span>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include "base_agent_container.h"
#include "base_agent_data.h"
#include "generic_agent_solver.h"

using namespace physicore;

TEST(BaseAgentContainerTest, CreateAndRemove)
{
	base_agent_container container(std::make_unique<base_agent_data>());
	auto* agent1 = container.create();
	EXPECT_EQ(container.size(), 1);
	EXPECT_NE(agent1, nullptr);
	auto* agent2 = container.create();
	EXPECT_EQ(container.size(), 2);
	EXPECT_NE(agent2, nullptr);
	EXPECT_NE(agent1, agent2);
	container.remove_agent(agent1);
	EXPECT_EQ(container.size(), 1);
	// After removal, agent1 pointer is invalid, but agent2 should still be valid
	EXPECT_NE(container.create(), nullptr);
}

class RemoveAgentTest : public ::testing::TestWithParam<int>
{};

TEST_P(RemoveAgentTest, RemoveAgentsAndCheckPositions)
{
	base_agent_container container(std::make_unique<base_agent_data>());
	auto* agent0 = container.create();
	auto* agent1 = container.create();
	auto* agent2 = container.create();

	agent0->position()[0] = 1.0;
	agent0->position()[1] = 2.0; // agent 0
	agent1->position()[0] = 3.0;
	agent1->position()[1] = 4.0; // agent 1
	agent2->position()[0] = 5.0;
	agent2->position()[1] = 6.0; // agent 2

	const int remove_idx = GetParam();
	container.remove_at(remove_idx);

	if (remove_idx != 0)
	{
		EXPECT_EQ(agent0->position()[0], 1.0);
		EXPECT_EQ(agent0->position()[1], 2.0);
	}
	if (remove_idx != 1)
	{
		EXPECT_EQ(agent1->position()[0], 3.0);
		EXPECT_EQ(agent1->position()[1], 4.0);
	}
	if (remove_idx != 2)
	{
		EXPECT_EQ(agent2->position()[0], 5.0);
		EXPECT_EQ(agent2->position()[1], 6.0);
	}
}

INSTANTIATE_TEST_SUITE_P(BaseAgentContainerTest, RemoveAgentTest, ::testing::Values(0, 1, 2));

class diffusion_agent_data
{
public:
	base_agent_data& data;

	void add() {}

	void remove_at(index_t /* position */) {}

	diffusion_agent_data(base_agent_data& data) : data(data) {}
};

class diffusion_agent_interface : public virtual base_agent_interface
{};

class diffusion_agent : public base_agent, public virtual diffusion_agent_interface
{
	diffusion_agent_data& data;

public:
	using DataType = diffusion_agent_data;
	using InterfaceType = diffusion_agent_interface;

	diffusion_agent(index_t index,
					std::tuple<std::unique_ptr<base_agent_data>, std::unique_ptr<diffusion_agent_data>>& datas)
		: base_agent_interface(index), base_agent(index, *std::get<0>(datas)), data(*std::get<1>(datas))
	{}
};

class mechanics_agent_data
{
public:
	base_agent_data& data;

	void add() {}

	void remove_at(index_t /* position */) {}

	mechanics_agent_data(base_agent_data& data) : data(data) {}
};

class mechanics_agent_interface : public virtual base_agent_interface
{};

class mechanics_agent : public base_agent, public virtual mechanics_agent_interface
{
	mechanics_agent_data& data;

public:
	using DataType = mechanics_agent_data;
	using InterfaceType = mechanics_agent_interface;

	mechanics_agent(index_t index,
					std::tuple<std::unique_ptr<base_agent_data>, std::unique_ptr<mechanics_agent_data>>& datas)
		: base_agent_interface(index), base_agent(index, *std::get<0>(datas)), data(*std::get<1>(datas))
	{}
};

class big_agent_data
{
public:
	base_agent_data& data1;
	diffusion_agent_data& data2;
	mechanics_agent_data& data3;

	void add() {}

	void remove_at(index_t /* position */) {}

	big_agent_data(base_agent_data& data1, diffusion_agent_data& data2, mechanics_agent_data& data3)
		: data1(data1), data2(data2), data3(data3)
	{}
};

class big_agent_interface : public virtual base_agent_interface,
							public virtual diffusion_agent_interface,
							public virtual mechanics_agent_interface
{};

class big_agent : public base_agent, public virtual big_agent_interface
{
public:
	big_agent_data& data;

	using DataType = big_agent_data;
	using InterfaceType = big_agent_interface;

	big_agent(index_t index, std::tuple<std::unique_ptr<base_agent_data>, std::unique_ptr<diffusion_agent_data>,
										std::unique_ptr<mechanics_agent_data>, std::unique_ptr<big_agent_data>>& datas)
		: base_agent_interface(index), base_agent(index, *std::get<0>(datas)), data(*std::get<3>(datas))
	{}
};

class diff_retriever : public generic_agent_solver<diffusion_agent>
{
public:
	using generic_agent_solver<diffusion_agent>::retrieve_agent_data;
};

class mech_retriever : public generic_agent_solver<mechanics_agent>
{
public:
	using generic_agent_solver<mechanics_agent>::retrieve_agent_data;
};

// An agent that is not built on base_agent: it has its own data and its own interface, as a module
// may need for objects that are not point agents of the shared kind.
struct standalone_data
{
	std::vector<real_t> positions;
	std::vector<real_t> values;

	void add()
	{
		positions.resize(positions.size() + 3, 0.0);
		values.push_back(0.0);
	}

	void remove_at(index_t position)
	{
		const index_t last = values.size() - 1;
		if (position != last)
		{
			values[position] = values[last];
			for (index_t d = 0; d < 3; ++d)
			{
				positions[3 * position + d] = positions[3 * last + d];
			}
		}
		values.pop_back();
		positions.resize(positions.size() - 3);
	}
};

class standalone_interface : public virtual base_agent_interface
{
public:
	virtual real_t& value() = 0;
};

class standalone_agent : public virtual standalone_interface
{
	standalone_data& data;

public:
	using DataType = standalone_data;
	using InterfaceType = standalone_interface;

	standalone_agent(index_t index, std::tuple<std::unique_ptr<standalone_data>>& datas)
		: base_agent_interface(index), data(*std::get<0>(datas))
	{}

	std::span<real_t> position() override { return { &data.positions[3 * this->index], 3 }; }

	real_t& value() override { return data.values[this->index]; }
};

TEST(BaseAgentContainerTest, ContainerOfAgentsNotBuiltOnBaseAgentRemovesAndKeepsIndicesValid)
{
	// The container hands out its agents through their own interface, here standalone_interface, and
	// base_agent is not among its types. Removing an agent moves the last one into its place, for which
	// the container must reach the index of an agent seen through that interface.
	generic_agent_and_data_container<standalone_agent> container(std::make_unique<standalone_data>());

	standalone_interface* first = container.create();
	standalone_interface* second = container.create();
	standalone_interface* third = container.create();
	first->value() = 1.0;
	second->value() = 2.0;
	third->value() = 3.0;
	third->position()[1] = 30.0;

	container.remove_agent(second);

	ASSERT_EQ(container.size(), 2U);
	EXPECT_EQ(container.get_agent_at(0), first);
	// the last agent took over the freed index, and its data moved with it
	EXPECT_EQ(container.get_agent_at(1), third);
	EXPECT_DOUBLE_EQ(container.get_agent_at(0)->value(), 1.0);
	EXPECT_DOUBLE_EQ(container.get_agent_at(1)->value(), 3.0);
	EXPECT_DOUBLE_EQ(container.get_agent_at(1)->position()[1], 30.0);

	container.remove_at(0);

	ASSERT_EQ(container.size(), 1U);
	EXPECT_EQ(container.get_agent_at(0), third);
	EXPECT_DOUBLE_EQ(third->value(), 3.0);
	EXPECT_DOUBLE_EQ(third->position()[1], 30.0);
}

TEST(BaseAgentContainerTest, Instantiation)
{
	const generic_agent_and_data_container<base_agent> base_container(std::make_unique<base_agent_data>());

	// Instantiate diffusion container
	{
		auto base_data = std::make_unique<physicore::base_agent_data>();
		auto diffusion_data = std::make_unique<diffusion_agent_data>(*base_data);

		const generic_agent_and_data_container<base_agent, diffusion_agent> container(std::move(base_data),
																					  std::move(diffusion_data));
	}

	// Instantiate mechanics container
	{
		auto base_data = std::make_unique<physicore::base_agent_data>();
		auto mechanics_data = std::make_unique<mechanics_agent_data>(*base_data);

		const generic_agent_and_data_container<base_agent, mechanics_agent> container(std::move(base_data),
																					  std::move(mechanics_data));
	}

	// Instantiate big (union of diffusion and mechanics) container
	{
		auto base_data = std::make_unique<physicore::base_agent_data>();
		auto diffusion_data = std::make_unique<diffusion_agent_data>(*base_data);
		auto mechanics_data = std::make_unique<mechanics_agent_data>(*base_data);
		auto big_data = std::make_unique<big_agent_data>(*base_data, *diffusion_data, *mechanics_data);

		generic_agent_and_data_container<base_agent, diffusion_agent, mechanics_agent, big_agent> container(
			std::move(base_data), std::move(diffusion_data), std::move(mechanics_data), std::move(big_data));

		// Can be down-casted
		generic_agent_interface_container<diffusion_agent_interface>& diffusion_container = container;
		generic_agent_interface_container<mechanics_agent_interface>& mechanics_container = container;

		// Data can be accessed
		{
			diff_retriever().retrieve_agent_data(container);

			mech_retriever().retrieve_agent_data(container);
		}

		mechanics_container.create();
		ASSERT_EQ(diffusion_container.size(), 1);
		ASSERT_EQ(mechanics_container.size(), 1);

		diffusion_container.create();
		ASSERT_EQ(diffusion_container.size(), 2);
		ASSERT_EQ(mechanics_container.size(), 2);
	}
}
