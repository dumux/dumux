add_input_file_links()

dune_add_test(NAME test_pipe_laufer
              SOURCES test_pipe_laufer.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/pipe_laufer.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/pipe_laufer-00002.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_pipe_laufer")
