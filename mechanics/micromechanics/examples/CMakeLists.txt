set(ALL_EXAMPLES potential_curve two_cells)

foreach(example ${ALL_EXAMPLES})
  add_executable(mechanics.micromechanics.${example} ${example}.cpp)

  target_link_libraries(mechanics.micromechanics.${example}
                        PRIVATE mechanics.micromechanics common)
endforeach()
