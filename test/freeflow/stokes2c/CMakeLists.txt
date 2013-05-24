# add build targets
add_executable("test_stokes2c" test_stokes2c.cc)
add_dune_superlu_flags(test_stokes2c)
add_test(test_stokes
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/stokes2c-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/stokes2c-00007.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_stokes2c
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_stokes2c.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/test_stokes2c.dgf)

set(NORMALTESTS test_stokes2c)
set(ALLTESTS ${NORMALTESTS})

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${ALLTESTS})

foreach(_exe ${ALLTESTS})
  target_link_libraries(${_exe} ${DUNE_LIBS})
endforeach(_exe ${ALLTESTS})
