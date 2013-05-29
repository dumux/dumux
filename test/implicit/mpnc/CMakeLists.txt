# build target for the twophase-twocomponent test problem
add_executable("test_boxmpnc" EXCLUDE_FROM_ALL test_boxmpnc.cc)
add_test(test_boxmpnc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/obstaclebox-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/obstaclebox-00010.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_boxmpnc
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_boxmpnc.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/obstacle_24x16.dgf
  -TimeManager.TEnd 1e4
  -TimeManager.DtInitial 250)

add_executable("test_ccmpnc" EXCLUDE_FROM_ALL test_ccmpnc.cc)
add_test(test_ccmpnc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/obstaclecc-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/obstaclecc-00010.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_ccmpnc
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_ccmpnc.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/obstacle_24x16.dgf
  -TimeManager.TEnd 1e4
  -TimeManager.DtInitial 250)

# build target for the one-phase forchheimer test problem
add_executable("test_forchheimer1p" EXCLUDE_FROM_ALL test_forchheimer1p.cc)
add_test(test_forchheimer1p
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/forchheimer1p-reference.vtp
  ${CMAKE_CURRENT_BINARY_DIR}/test_forchheimer1p-00002.vtp
  ${CMAKE_CURRENT_BINARY_DIR}/test_forchheimer1p
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_forchheimer1p.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/forchheimer1d.dgf)

# build target for the two-phase forchheimer test problem
add_executable("test_forchheimer2p" EXCLUDE_FROM_ALL test_forchheimer2p.cc)
add_test(test_forchheimer2p
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/forchheimer2p-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_forchheimer2p-00009.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_forchheimer2p
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_forchheimer2p.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/obstacle_24x16.dgf)

set(NORMALTESTS
  test_boxmpnc
  test_ccmpnc
  test_forchheimer1p
  test_forchheimer2p)
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
