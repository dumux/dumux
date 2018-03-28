add_input_file_links()

dune_add_test(NAME test_channel_zeroeq2c
              SOURCES test_channel_zeroeq.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_channel_zeroeq2c.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_channel_zeroeq2c-00033.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_channel_zeroeq2c test_channel_zeroeq2c.input")
