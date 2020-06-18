dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME example_shallowwaterfriction
               LABELS shallowwater example
               SOURCES main.cc
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS       --script fuzzy
                              --files ${CMAKE_SOURCE_DIR}/test/references/example_shallowwaterfriction-reference.vtu
                                      ${CMAKE_CURRENT_BINARY_DIR}/example_shallowwaterfriction-00001.vtu
                              --zeroThreshold {"velocityY":1e-14}
                              --command "${CMAKE_CURRENT_BINARY_DIR}/example_shallowwaterfriction params.input
                              -Problem.Name example_shallowwaterfriction")
