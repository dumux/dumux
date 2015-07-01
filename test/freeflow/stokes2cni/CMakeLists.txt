add_dumux_test(test_stokes2cni test_stokes2cni test_stokes2cni.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/stokes2cni-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/stokes2cni-00008.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_stokes2cni)

add_dumux_test(test_stokes2cni_conduction test_stokes2cni test_stokes2cni.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/stokes2cni_conduction-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/stokes2cni_conduction-00008.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_stokes2cni
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_stokes2cni_conduction.input)
