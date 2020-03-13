dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME test_shallowwater_roughchannel
               SOURCES main.cc
               LABELS shallowwater
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS      --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_shallowwater_roughchannel-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/roughchannel-00001.vtu
                             --zeroThreshold {"velocityY":1e-14}
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_shallowwater_roughchannel params.input
                             -Problem.Name roughchannel")
