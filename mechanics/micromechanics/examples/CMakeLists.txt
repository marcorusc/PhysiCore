set(ALL_EXAMPLES
    cell_sorting
    domain_walls
    growth_division
    motility
    potential_curve
    pressure
    run_settings
    scaled_cell
    thermal_noise
    two_cells)

foreach(example ${ALL_EXAMPLES})
  add_executable(mechanics.micromechanics.${example} ${example}.cpp)

  target_link_libraries(mechanics.micromechanics.${example}
                        PRIVATE mechanics.micromechanics common)
endforeach()
