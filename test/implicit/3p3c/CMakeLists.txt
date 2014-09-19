add_dumux_test(test_box3p3c test_box3p3c test_box3p3c.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/infiltrationbox-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/box3p3c-00007.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_box3p3c
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_box3p3c_reference.input)

add_dumux_test(test_cc3p3c test_cc3p3c test_cc3p3c.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/infiltrationcc-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/cc3p3c-00006.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_cc3p3c
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_cc3p3c_reference.input)
