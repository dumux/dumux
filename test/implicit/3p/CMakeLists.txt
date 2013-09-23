add_dumux_test(test_box3p test_box3p test_box3p.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/infiltration3pbox-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/infiltration3pbox-00009.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_box3p
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_box3p.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/infiltration3p.dgf)

add_dumux_test(test_cc3p test_cc3p test_cc3p.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/infiltration3pcc-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/infiltration3pcc-00009.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_cc3p
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_cc3p.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/infiltration3p.dgf)
