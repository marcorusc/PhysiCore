add_library(mechanics.micromechanics)
add_library(mechanics::micromechanics ALIAS mechanics.micromechanics)
add_library(mechanics.micromechanics_internal_iface INTERFACE)
add_library(mechanics.micromechanics_iface INTERFACE)

include(GenerateExportHeader)
generate_export_header(
  mechanics.micromechanics BASE_NAME micromechanics EXPORT_FILE_NAME
  ${CMAKE_CURRENT_BINARY_DIR}/include/micromechanics/micromechanics_export.h)

if(${BUILD_SHARED_LIBS})
  set_target_properties(
    mechanics.micromechanics
    PROPERTIES CXX_VISIBILITY_PRESET hidden
               VISIBILITY_INLINES_HIDDEN ON
               POSITION_INDEPENDENT_CODE ON
               INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()

file(GLOB PUBLIC_INC "include/micromechanics/*.h")
target_sources(
  mechanics.micromechanics
  PRIVATE src/agent.cpp
          src/agent_data.cpp
          src/cell.cpp
          src/cell_data.cpp
          src/cell_parameters.cpp
          src/cell_placement.cpp
          src/cell_type.cpp
          src/config_reader.cpp
          src/environment.cpp
          src/environment_builder.cpp
          src/solver_registry.cpp
          src/vtk_agents_serializer.cpp
          ${PROJECT_SOURCE_DIR}/common/src/vtk_serializer_base.cpp
  PUBLIC
    FILE_SET
    HEADERS
    BASE_DIRS
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "${CMAKE_CURRENT_BINARY_DIR}/include"
    FILES
    ${PUBLIC_INC}
    "${CMAKE_CURRENT_BINARY_DIR}/include/micromechanics/micromechanics_export.h"
)

target_include_directories(
  mechanics.micromechanics_internal_iface
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/micromechanics>)

# what the kernels compile against: the public headers, without import or export marks, as the kernels
# become part of the library itself
target_sources(
  mechanics.micromechanics_iface
  INTERFACE FILE_SET HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
            "${CMAKE_CURRENT_BINARY_DIR}/include")

target_compile_definitions(mechanics.micromechanics_iface
                           INTERFACE MICROMECHANICS_STATIC_DEFINE)

target_link_libraries(mechanics.micromechanics_internal_iface INTERFACE common)

target_link_libraries(mechanics.micromechanics_iface INTERFACE common)

find_package(pugixml CONFIG REQUIRED)
find_package(VTK CONFIG REQUIRED)

target_link_libraries(mechanics.micromechanics_internal_iface
                      INTERFACE pugixml::pugixml VTK::IOXML)

target_link_libraries(
  mechanics.micromechanics
  PRIVATE mechanics.micromechanics_internal_iface
  PUBLIC common)

add_subdirectory(kernels/openmp_solver)

if(PHYSICORE_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

if(PHYSICORE_BUILD_TESTS)
  add_subdirectory(tests)
endif()

install(
  TARGETS mechanics.micromechanics mechanics.micromechanics_internal_iface
  EXPORT MicromechanicsTargets
  FILE_SET HEADERS)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/MicromechanicsConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)

configure_package_config_file(
  cmake/MicromechanicsConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/MicromechanicsConfig.cmake"
  INSTALL_DESTINATION lib/cmake/Micromechanics)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/MicromechanicsConfig.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/MicromechanicsConfigVersion.cmake"
        DESTINATION lib/cmake/Micromechanics)

install(
  EXPORT MicromechanicsTargets
  NAMESPACE PhysiCore::
  FILE MicromechanicsTargets.cmake
  DESTINATION lib/cmake/Micromechanics)
