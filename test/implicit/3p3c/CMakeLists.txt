add_dumux_test(test_box3p3c test_box3p3c test_box3p3c.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/infiltrationbox-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/infiltrationbox-00007.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_box3p3c
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_box3p3c.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/test_3p3c_coarse.dgf
  -TimeManager.TEnd 8.64e5
  -TimeManager.DtInitial 8.64e4)

add_dumux_test(test_cc3p3c test_cc3p3c test_cc3p3c.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/infiltrationcc-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/infiltrationcc-00006.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_cc3p3c
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_cc3p3c.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/test_3p3c_coarse.dgf
  -TimeManager.TEnd 8.64e5
  -TimeManager.DtInitial 8.64e4)
