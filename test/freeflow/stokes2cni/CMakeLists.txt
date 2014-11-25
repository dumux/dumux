add_dumux_test(test_stokes2cni test_stokes2cni test_stokes2cni.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/stokes2cni-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/stokes2cni-00008.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_stokes2cni)
