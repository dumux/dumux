add_input_file_links()

add_dumux_test(test_2pdfm test_2pdfm test_2pdfm.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/2pdfm-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/2pdfm-00011.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pdfm")

#install sources
install(FILES
2pdfmspatialparams.hh
2pdfmtestproblem.hh
test_2pdfm.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/2pdfm)
