add_executable("test_box1p2c" test_box1p2c.cc)
add_test(test_box1p2c
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/outflowbox-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/outflowbox-00010.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_box1p2c
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_box1p2c.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/test_1p2c.dgf
  -TimeManager.TEnd 100
  -TimeManager.DtInitial 1)

add_executable("test_cc1p2c" test_cc1p2c.cc)
add_test(test_cc1p2c
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/outflowcc-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/outflowcc-00010.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_cc1p2c
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_cc1p2c.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/test_1p2c.dgf
  -TimeManager.TEnd 100
  -TimeManager.DtInitial 1)

set(NORMALTESTS
  test_box1p2c
  test_cc1p2c)
set(ALLTESTS ${NORMALTESTS})

# We do not want want to build the tests during make all,
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${ALLTESTS})

foreach(_exe ${ALLTESTS})
  target_link_libraries(${_exe} ${DUNE_LIBS})
endforeach(_exe ${ALLTESTS})
