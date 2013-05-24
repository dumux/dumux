# add build targets
add_executable("test_navierstokes" test_navierstokes.cc)
add_dune_superlu_flags(test_navierstokes)
add_dune_alugrid_flags(test_navierstokes)
add_dune_ug_flags(test_navierstokes)
add_test(test_navierstokes
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/navierstokes-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/navierstokes-00008.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_navierstokes
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_navierstokes.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/test_navierstokes.dgf)

set(NORMALTESTS test_navierstokes)
set(ALLTESTS ${NORMALTESTS})

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${ALLTESTS})

foreach(_exe ${ALLTESTS})
  target_link_libraries(${_exe} ${DUNE_LIBS})
endforeach(_exe ${ALLTESTS})
