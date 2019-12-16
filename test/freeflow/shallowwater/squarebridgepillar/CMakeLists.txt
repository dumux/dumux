dune_symlink_to_source_files(FILES "params.input")

dune_add_test(NAME test_shallowwater_squarebridgepillar
              SOURCES main.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_shallowwater_squarebridgepillar-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/squarebridgepillar-00001.vtu
                             --zeroThreshold {"velocityY":1e-14}
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_shallowwater_squarebridgepillar params.input
                             -Problem.Name squarebridgepillar")
