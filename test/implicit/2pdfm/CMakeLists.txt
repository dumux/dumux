add_dumux_test(test_2pdfm test_2pdfm test_2pdfm.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/2pdfm-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/2pdfm-00011.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pdfm")
