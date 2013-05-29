# build target for the 2p2cni test problem
add_executable("test_box2p2cni" EXCLUDE_FROM_ALL test_box2p2cni.cc)
add_test(test_box2p2cni
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/waterairbox-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/waterairbox-00010.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_box2p2cni
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_box2p2cni.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/test_2p2cni.dgf
  -TimeManager.TEnd 1e4
  -TimeManager.DtInitial 250)

add_executable("test_cc2p2cni" EXCLUDE_FROM_ALL test_cc2p2cni.cc)
add_test(test_cc2p2cni
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/wateraircc-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/wateraircc-00010.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_cc2p2cni
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_cc2p2cni.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/test_2p2cni.dgf
  -TimeManager.TEnd 1e4
  -TimeManager.DtInitial 250)

set(NORMALTESTS
  test_box2p2cni
  test_cc2p2cni)
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
