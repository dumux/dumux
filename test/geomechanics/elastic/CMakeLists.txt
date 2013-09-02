# add build targets
add_executable("test_elastic" EXCLUDE_FROM_ALL test_elastic.cc)
add_dune_superlu_flags(test_elastic)
add_test(test_elastic
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/elasticmatrix-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/elasticmatrix-00002.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_elastic
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_elastic.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/test_elastic.dgf)

set(NORMALTESTS test_elastic)
set(ALLTESTS ${NORMALTESTS})

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${ALLTESTS})
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  PROPERTY TEST_INCLUDE_FILE ${CMAKE_CURRENT_BINARY_DIR}/BuildTests.cmake)

foreach(_exe ${ALLTESTS})
  target_link_libraries(${_exe} ${DUNE_LIBS})
endforeach(_exe ${ALLTESTS})
