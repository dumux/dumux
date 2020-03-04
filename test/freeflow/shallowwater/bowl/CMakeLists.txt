dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME test_shallowwater_bowl
              SOURCES main.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_shallowwater_bowl-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/bowl-00001.vtu
                             --zeroThreshold {"velocityY":1e-14}
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_shallowwater_bowl params.input
                             -Problem.Name bowl")
