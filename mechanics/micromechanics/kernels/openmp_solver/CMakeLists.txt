add_library(mechanics.micromechanics.kernels.openmp_solver_internal_iface
            INTERFACE)
add_library(mechanics.micromechanics.kernels.openmp_solver OBJECT)

if(${BUILD_SHARED_LIBS})
  set_target_properties(mechanics.micromechanics.kernels.openmp_solver
                        PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

find_package(OpenMP 4)

target_sources(
  mechanics.micromechanics.kernels.openmp_solver
  PRIVATE src/cell_solver.cpp
          src/force_solver.cpp
          src/neighbor_grid.cpp
          src/openmp_solver.cpp
          src/position_solver.cpp
          src/register_solver.cpp
  PUBLIC FILE_SET HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_link_libraries(
  mechanics.micromechanics.kernels.openmp_solver_internal_iface
  INTERFACE $<BUILD_INTERFACE:mechanics.micromechanics_iface>)

target_include_directories(
  mechanics.micromechanics.kernels.openmp_solver_internal_iface
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>)

if(OpenMP_CXX_FOUND)
  target_link_libraries(
    mechanics.micromechanics.kernels.openmp_solver_internal_iface
    INTERFACE OpenMP::OpenMP_CXX)
endif()

target_link_libraries(
  mechanics.micromechanics.kernels.openmp_solver
  PRIVATE mechanics.micromechanics.kernels.openmp_solver_internal_iface)

target_link_libraries(mechanics.micromechanics
                      PRIVATE mechanics.micromechanics.kernels.openmp_solver)

if(PHYSICORE_BUILD_TESTS)
  add_subdirectory(tests)
endif()
