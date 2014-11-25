add_dumux_test(test_el1p2c test_el1p2c test_el1p2c.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/el1p2c-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/el1p2c-00015.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_el1p2c)
