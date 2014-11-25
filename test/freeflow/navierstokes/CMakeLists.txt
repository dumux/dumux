add_dumux_test(test_navierstokes test_navierstokes test_navierstokes.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/navierstokes-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/navierstokes-00008.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_navierstokes)
