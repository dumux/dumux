add_input_file_links()
dune_symlink_to_source_files(FILES laufer_re50000_u+y+.csv)

dune_add_test(NAME test_pipe_laufer
              SOURCES test_pipe_laufer.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/pipe_laufer_zeroeq.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/pipe_laufer_zeroeq_reference-00021.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_pipe_laufer test_pipe_laufer_reference.input")
