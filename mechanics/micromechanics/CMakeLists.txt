add_library(mechanics.micromechanics STATIC)
add_library(mechanics::micromechanics ALIAS mechanics.micromechanics)
add_library(mechanics.micromechanics_internal_iface INTERFACE)

file(GLOB PUBLIC_INC "include/micromechanics/*.h")
target_sources(
  mechanics.micromechanics
  PRIVATE src/agent.cpp
          src/agent_data.cpp
          src/cell.cpp
          src/cell_data.cpp
          src/environment.cpp
          src/solver_registry.cpp
          src/vtk_agents_serializer.cpp
          ${PROJECT_SOURCE_DIR}/common/src/vtk_serializer_base.cpp
  PUBLIC FILE_SET HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include" FILES
         ${PUBLIC_INC})

target_sources(
  mechanics.micromechanics_internal_iface
  INTERFACE FILE_SET HEADERS BASE_DIRS
            "${CMAKE_CURRENT_SOURCE_DIR}/include/micromechanics"
            "${CMAKE_CURRENT_SOURCE_DIR}/src")

target_link_libraries(mechanics.micromechanics_internal_iface INTERFACE common)

find_package(VTK CONFIG REQUIRED)

target_link_libraries(mechanics.micromechanics_internal_iface
                      INTERFACE VTK::IOXML)

add_library(mechanics.micromechanics_iface INTERFACE)
target_include_directories(
  mechanics.micromechanics_iface
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(mechanics.micromechanics_iface INTERFACE common)

target_link_libraries(mechanics.micromechanics
                      PRIVATE mechanics.micromechanics_internal_iface)

add_subdirectory(kernels/openmp_solver)

if(PHYSICORE_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

if(PHYSICORE_BUILD_TESTS)
  add_subdirectory(tests)
endif()
