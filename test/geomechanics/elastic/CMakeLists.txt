add_input_file_links()

add_dumux_test(test_elastic test_elastic test_elastic.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/elasticmatrix-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/elasticmatrix-00002.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_elastic")
