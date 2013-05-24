# add build targets
add_executable("tutorial_decoupled" tutorial_decoupled.cc)
add_test(tutorial_decoupled
  ${CMAKE_CURRENT_BINARY_DIR}/tutorial_decoupled
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/tutorial_decoupled.input)

add_executable("tutorial_coupled" tutorial_coupled.cc)
add_dune_alugrid_flags(tutorial_coupled)
add_test(tutorial_coupled
  ${CMAKE_CURRENT_BINARY_DIR}/tutorial_coupled
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/tutorial_coupled.input)

set(NORMALTESTS
  tutorial_decoupled
  tutorial_coupled)
set(ALLTESTS ${NORMALTESTS})

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${ALLTESTS})

foreach(_exe ${ALLTESTS})
  target_link_libraries(${_exe} ${DUNE_LIBS})
endforeach(_exe ${ALLTESTS})
