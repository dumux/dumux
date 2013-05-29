# build target for the simple twophase lens problem
add_executable("test_box2p" EXCLUDE_FROM_ALL test_box2p.cc)
add_test(test_box2p
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/lensbox-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/lensbox-00009.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_box2p
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_box2p.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/test_2p.dgf
  -TimeManager.TEnd 3000
  -TimeManager.DtInitial 250)

add_executable("test_cc2p" EXCLUDE_FROM_ALL test_cc2p.cc)
add_test(test_cc2p
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/lenscc-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/lenscc-00009.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_cc2p
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_cc2p.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/test_2p.dgf
  -TimeManager.TEnd 3000
  -TimeManager.DtInitial 250)

set(NORMALTESTS
  test_box2p
  test_cc2p)
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
