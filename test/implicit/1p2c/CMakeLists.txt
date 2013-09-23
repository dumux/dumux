add_dumux_test(test_box1p2c test_box1p2c test_box1p2c.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/outflowbox-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/outflowbox-00010.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_box1p2c
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_box1p2c.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/test_1p2c.dgf
  -TimeManager.TEnd 100
  -TimeManager.DtInitial 1)

add_dumux_test(test_cc1p2c test_cc1p2c test_cc1p2c.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/outflowcc-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/outflowcc-00010.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_cc1p2c
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_cc1p2c.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/test_1p2c.dgf
  -TimeManager.TEnd 100
  -TimeManager.DtInitial 1)
