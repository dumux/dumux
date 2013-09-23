add_dumux_test(test_2pdfm test_2pdfm test_2pdfm.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/2pdfm-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/2pdfm-00011.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_2pdfm
  -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_2pdfm.input
  -Grid.File ${CMAKE_CURRENT_SOURCE_DIR}/grids/2pdfmartmesh.net)
add_dune_ug_flags(test_2pdfm)
add_dune_mpi_flags(test_2pdfm)
