#include <gtest/gtest.h>

#include "micromechanics/solver_registry.h"

namespace physicore::mechanics::micromechanics::tests {

namespace {
class registered_test_solver : public solver
{
public:
	void initialize(environment& /* e */) override {}
	void solve(environment& /* e */, index_t /* iterations */) override {}
};
} // namespace

TEST(SolverRegistryTest, UnknownSolverReturnsNullptr)
{
	EXPECT_EQ(solver_registry::instance().get("nonexistent_solver"), nullptr);
}

TEST(SolverRegistryTest, RegisteredSolverCanBeInstantiated)
{
	const registry_adder<registered_test_solver> adder("micromechanics_test_solver");

	auto instance = solver_registry::instance().get("micromechanics_test_solver");

	ASSERT_NE(instance, nullptr);
	EXPECT_NE(dynamic_cast<registered_test_solver*>(instance.get()), nullptr);
}

} // namespace physicore::mechanics::micromechanics::tests
