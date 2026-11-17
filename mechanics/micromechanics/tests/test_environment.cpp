#include <memory>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "micromechanics/environment.h"

namespace physicore::mechanics::micromechanics::tests {

namespace {
class mock_solver : public solver
{
public:
	void initialize(environment& /* e */) override {}
	void solve(environment& /* e */, index_t iterations) override { solved_iterations += iterations; }
	index_t solved_iterations = 0;
};

class mock_serializer : public serializer
{
public:
	void serialize(const environment& /* e */, real_t current_time) override { serialized_time = current_time; }
	real_t serialized_time = -1.0;
};

environment make_environment() { return { { 3, { -100, -100, -100 }, { 100, 100, 100 }, { 20, 20, 20 } }, 0.1 }; }

// Cell index of every sub-agent, in container order
std::vector<index_t> agent_cell_indices(environment& env)
{
	std::vector<index_t> indices;
	for (index_t i = 0; i < env.agents->size(); ++i)
	{
		indices.push_back(env.agents->get_agent_at(i)->cell_index());
	}
	return indices;
}

// Number of sub-agents of every cell, in container order
std::vector<index_t> cell_agent_counts(environment& env)
{
	std::vector<index_t> counts;
	for (index_t i = 0; i < env.cells->size(); ++i)
	{
		counts.push_back(env.cells->get_agent_at(i)->agent_count());
	}
	return counts;
}
} // namespace

TEST(EnvironmentTest, ConstructorCreatesEmptyContainers)
{
	const environment env = make_environment();

	EXPECT_DOUBLE_EQ(env.mechanics_timestep, 0.1);
	EXPECT_EQ(env.mesh.dims, index_t { 3 });
	ASSERT_NE(env.cells, nullptr);
	ASSERT_NE(env.agents, nullptr);
	EXPECT_EQ(env.cells->size(), 0U);
	EXPECT_EQ(env.agents->size(), 0U);
	EXPECT_EQ(env.solver, nullptr);
	EXPECT_EQ(env.serializer, nullptr);
}

TEST(EnvironmentTest, RunSingleTimestepAndSerializeStateDelegate)
{
	environment env = make_environment();

	EXPECT_THROW(env.run_single_timestep(), std::runtime_error);
	EXPECT_NO_THROW(env.serialize_state(1.0));

	auto solver_owner = std::make_unique<mock_solver>();
	auto serializer_owner = std::make_unique<mock_serializer>();
	const auto* used_solver = solver_owner.get();
	const auto* used_serializer = serializer_owner.get();
	env.solver = std::move(solver_owner);
	env.serializer = std::move(serializer_owner);

	env.run_single_timestep();
	env.run_single_timestep();
	EXPECT_EQ(used_solver->solved_iterations, index_t { 2 });

	env.serialize_state(1.5);
	EXPECT_DOUBLE_EQ(used_serializer->serialized_time, 1.5);
}

TEST(EnvironmentTest, CreateAgentAssignsCellAndCountsIt)
{
	environment env = make_environment();

	auto* cell0 = env.create_cell();
	auto* cell1 = env.create_cell();
	ASSERT_EQ(env.cells->size(), 2U);
	EXPECT_EQ(cell0->agent_count(), index_t { 0 });

	auto* a0 = env.create_agent(1);
	auto* a1 = env.create_agent(0);
	auto* a2 = env.create_agent(1);

	EXPECT_EQ(a0->cell_index(), index_t { 1 });
	EXPECT_EQ(a1->cell_index(), index_t { 0 });
	EXPECT_EQ(a2->cell_index(), index_t { 1 });
	EXPECT_EQ(a0->position().size(), 3U);
	EXPECT_EQ(cell0->agent_count(), index_t { 1 });
	EXPECT_EQ(cell1->agent_count(), index_t { 2 });

	EXPECT_THROW(env.create_agent(2), std::out_of_range);
	EXPECT_EQ(env.agents->size(), 3U);
}

TEST(EnvironmentTest, RemoveAgentUpdatesCellCount)
{
	environment env = make_environment();
	env.create_cell();
	env.create_cell();
	env.create_agent(0);
	env.create_agent(1);
	auto* last = env.create_agent(1);
	last->position()[0] = 9.0;

	env.remove_agent(0);

	EXPECT_EQ(cell_agent_counts(env), std::vector<index_t>({ 0, 2 }));
	EXPECT_EQ(agent_cell_indices(env), std::vector<index_t>({ 1, 1 }));
	// the last sub-agent took over the freed index
	EXPECT_EQ(env.agents->get_agent_at(0), last);
	EXPECT_DOUBLE_EQ(last->position()[0], 9.0);

	EXPECT_THROW(env.remove_agent(2), std::out_of_range);
}

TEST(EnvironmentTest, RemoveCellRemovesItsAgentsAndRelinksMovedCell)
{
	environment env = make_environment();

	// three cells with 2, 1 and 3 sub-agents, created interleaved
	for (index_t i = 0; i < 3; ++i)
	{
		env.create_cell()->agent_radius() = static_cast<real_t>(i + 1);
	}
	for (const index_t cell_index : { 0, 2, 1, 2, 0, 2 })
	{
		env.create_agent(cell_index)->velocity()[0] = static_cast<real_t>(10 * (cell_index + 1));
	}

	env.remove_cell(0);

	// cell 2 (radius 3) moved to index 0; cell 1 stays
	ASSERT_EQ(env.cells->size(), 2U);
	EXPECT_DOUBLE_EQ(env.cells->get_agent_at(0)->agent_radius(), 3.0);
	EXPECT_DOUBLE_EQ(env.cells->get_agent_at(1)->agent_radius(), 2.0);
	EXPECT_EQ(cell_agent_counts(env), std::vector<index_t>({ 3, 1 }));

	// only the sub-agents of the removed cell are gone and every survivor points at its own cell
	ASSERT_EQ(env.agents->size(), 4U);
	index_t agents_of_moved_cell = 0;
	for (index_t i = 0; i < env.agents->size(); ++i)
	{
		auto* a = env.agents->get_agent_at(i);
		if (a->cell_index() == 0)
		{
			EXPECT_DOUBLE_EQ(a->velocity()[0], 30.0);
			++agents_of_moved_cell;
		}
		else
		{
			EXPECT_EQ(a->cell_index(), index_t { 1 });
			EXPECT_DOUBLE_EQ(a->velocity()[0], 20.0);
		}
	}
	EXPECT_EQ(agents_of_moved_cell, index_t { 3 });
}

TEST(EnvironmentTest, RemoveLastCellLeavesOtherLinksUntouched)
{
	environment env = make_environment();
	env.create_cell();
	env.create_cell();
	env.create_agent(0);
	env.create_agent(1);
	env.create_agent(0);

	env.remove_cell(1);

	EXPECT_EQ(cell_agent_counts(env), std::vector<index_t>({ 2 }));
	EXPECT_EQ(agent_cell_indices(env), std::vector<index_t>({ 0, 0 }));

	EXPECT_THROW(env.remove_cell(1), std::out_of_range);

	env.remove_cell(0);
	EXPECT_EQ(env.cells->size(), 0U);
	EXPECT_EQ(env.agents->size(), 0U);
}

} // namespace physicore::mechanics::micromechanics::tests
