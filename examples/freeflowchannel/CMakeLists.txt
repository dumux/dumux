dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME example_freeflow_channel_navierstokes
               LABELS freeflow navierstokes example
               SOURCES main.cc
               CMAKE_GUARD HAVE_UMFPACK
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS       --script fuzzy
                              --files ${CMAKE_SOURCE_DIR}/test/references/example_ff_channel-reference.vtu
                                      ${CMAKE_CURRENT_BINARY_DIR}/example_ff_channel-00001.vtu
                              --command "${CMAKE_CURRENT_BINARY_DIR}/example_freeflow_channel_navierstokes params.input -Problem.Name example_ff_channel")
