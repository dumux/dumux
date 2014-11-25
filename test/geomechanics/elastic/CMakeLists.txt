add_dumux_test(test_elastic test_elastic test_elastic.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/elasticmatrix-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/elasticmatrix-00002.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_elastic)
