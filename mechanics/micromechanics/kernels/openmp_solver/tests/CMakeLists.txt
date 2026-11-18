find_package(GTest REQUIRED)

file(GLOB TEST_SOURCES *.cpp)
add_executable(mechanics.micromechanics.kernels.openmp_solver.tests
               ${TEST_SOURCES})

target_link_libraries(
  mechanics.micromechanics.kernels.openmp_solver.tests
  PRIVATE mechanics.micromechanics.kernels.openmp_solver_internal_iface
          mechanics.micromechanics GTest::gtest GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(mechanics.micromechanics.kernels.openmp_solver.tests)
