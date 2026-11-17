set(ALL_EXAMPLES
    domain_walls
    growth_division
    motility
    potential_curve
    pressure
    scaled_cell
    thermal_noise
    two_cells)

foreach(example ${ALL_EXAMPLES})
  add_executable(mechanics.micromechanics.${example} ${example}.cpp)

  target_link_libraries(mechanics.micromechanics.${example}
                        PRIVATE mechanics.micromechanics common)
endforeach()
