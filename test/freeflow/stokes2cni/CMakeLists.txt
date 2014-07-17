add_dumux_test(test_stokes2cni test_stokes2cni test_stokes2cni.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/stokes2cni-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/stokes2cni-00008.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_stokes2cni
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_stokes2cni.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/test_stokes2cni.dgf)
add_dune_superlu_flags(test_stokes2cni)
add_dune_umfpack_flags(test_stokes2cni)
