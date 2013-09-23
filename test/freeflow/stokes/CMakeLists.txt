add_dumux_test(test_stokes test_stokes test_stokes.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/stokes-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/stokes-00014.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_stokes
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_stokes.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/test_stokes.dgf)
add_dune_superlu_flags(test_stokes)
