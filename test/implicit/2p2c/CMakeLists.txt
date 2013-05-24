# build target for the twophase-twocomponent test problem
add_executable("test_box2p2c" test_box2p2c.cc)
add_test(test_box2p2c
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/injectionbox-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/injectionbox-00009.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_box2p2c
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_box2p2c.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/test_2p2c.dgf
  -TimeManager.TEnd 1e4
  -TimeManager.DtInitial 250)

add_executable("test_cc2p2c" test_cc2p2c.cc)
add_test(test_cc2p2c
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/injectioncc-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/injectioncc-00009.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_cc2p2c
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_cc2p2c.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/test_2p2c.dgf
  -TimeManager.TEnd 1e4
  -TimeManager.DtInitial 250)

set(NORMALTESTS
  test_box2p2c
  test_cc2p2c)
set(ALLTESTS ${NORMALTESTS})

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${ALLTESTS})

foreach(_exe ${ALLTESTS})
  target_link_libraries(${_exe} ${DUNE_LIBS})
endforeach(_exe ${ALLTESTS})
