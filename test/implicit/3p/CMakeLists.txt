add_executable("test_box3p" EXCLUDE_FROM_ALL test_box3p.cc)
add_test(test_box3p
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/infiltration3pbox-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/infiltration3pbox-00009.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_box3p
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_box3p.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/infiltration3p.dgf)

add_executable("test_cc3p" EXCLUDE_FROM_ALL test_cc3p.cc)
add_test(test_cc3p
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/infiltration3pcc-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/infiltration3pcc-00009.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_cc3p
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_cc3p.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/infiltration3p.dgf)


set(NORMALTESTS
  test_box3p
  test_cc3p)
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
