add_dumux_test(test_stokes test_stokes test_stokes.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/stokes-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/stokes-00014.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_stokes)
