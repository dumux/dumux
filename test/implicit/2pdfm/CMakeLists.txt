add_executable("test_2pdfm" EXCLUDE_FROM_ALL test_2pdfm.cc)
add_dune_ug_flags(test_2pdfm)
add_test(test_2pdfm
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/2pdfm-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/2pdfm-00011.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_2pdfm
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_2pdfm.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/2pdfmartmesh.net)

set(NORMALTESTS test_2pdfm)
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
