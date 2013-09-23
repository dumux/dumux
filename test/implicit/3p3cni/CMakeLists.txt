add_dumux_test(test_box3p3cni test_box3p3cni test_box3p3cni.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/column3p3cnibox-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/columnxylolbox-00060.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_box3p3cni
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_box3p3cni.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/column.dgf)

add_dumux_test(test_cc3p3cni test_cc3p3cni test_cc3p3cni.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/column3p3cnicc-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/columnxylolcc-00054.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_cc3p3cni
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_cc3p3cni.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/column.dgf)
