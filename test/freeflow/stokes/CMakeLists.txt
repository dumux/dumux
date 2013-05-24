# add build targets
add_executable("test_stokes" test_stokes.cc)
add_dune_superlu_flags(test_stokes)
add_test(test_stokes
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/stokes-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/stokes-00013.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_stokes
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_stokes.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/test_stokes.dgf)

set(NORMALTESTS test_stokes)
set(ALLTESTS ${NORMALTESTS})

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${ALLTESTS})

foreach(_exe ${ALLTESTS})
  target_link_libraries(${_exe} ${DUNE_LIBS})
endforeach(_exe ${ALLTESTS})
