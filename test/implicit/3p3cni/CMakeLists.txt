add_executable("test_box3p3cni" EXCLUDE_FROM_ALL test_box3p3cni.cc)
add_test(test_box3p3cni
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/column3p3cnibox-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/columnxylolbox-00060.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_box3p3cni
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_box3p3cni.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/column.dgf)

add_executable("test_cc3p3cni" EXCLUDE_FROM_ALL test_cc3p3cni.cc)
add_test(test_cc3p3cni
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/column3p3cnicc-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/columnxylolcc-00054.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_cc3p3cni
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_cc3p3cni.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/column.dgf)

set(NORMALTESTS
  test_box3p3cni
  test_cc3p3cni)
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
