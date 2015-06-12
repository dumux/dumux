# test: test_zeroeq2cni
add_dumux_test(test_zeroeq2cni test_zeroeq2cni test_zeroeq2cni.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/zeroeq2cni-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/zeroeq2cni-00011.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_zeroeq2cni
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_zeroeq2cni.input
  )
