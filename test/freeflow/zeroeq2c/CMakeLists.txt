# test: test_zeroeq2c
add_dumux_test(test_zeroeq2c test_zeroeq2c test_zeroeq2c.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/zeroeq2c-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/zeroeq2c-00013.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_zeroeq2c
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_zeroeq2c.input
  )
