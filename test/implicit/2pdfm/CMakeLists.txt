add_dumux_test(test_2pdfm test_2pdfm test_2pdfm.cc
  ${CMAKE_SOURCE_DIR}/bin/runTest.sh
  ${CMAKE_SOURCE_DIR}/bin/fuzzycomparevtu.py
  ${CMAKE_SOURCE_DIR}/test/references/2pdfm-reference.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/2pdfm-00011.vtu
  ${CMAKE_CURRENT_BINARY_DIR}/test_2pdfm)
target_link_libraries(test_2pdfm ${DUNE_LIBS}  ${MPI_LIBRARIES})
add_dune_ug_flags(test_2pdfm)
