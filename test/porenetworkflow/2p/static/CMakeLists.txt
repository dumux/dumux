add_input_file_links()

dune_add_test(NAME test_pnm_2p_static
              SOURCES main.cc
              LABELS porenetworkflow
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzyData
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_pnm_2p_static-reference.txt
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p_static_pc-s-curve.txt
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pnm_2p_static")
