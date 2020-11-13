dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME test_shallowwater_bowl
               SOURCES main.cc
               LABELS shallowwater
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS   --script fuzzy
                          --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_shallowwater_bowl-reference.vtu
                                  ${CMAKE_CURRENT_BINARY_DIR}/bowl-00013.vtu
                          --zeroThreshold {"velocityY":1e-14}
                          --command "${CMAKE_CURRENT_BINARY_DIR}/test_shallowwater_bowl")

dumux_add_test(NAME test_shallowwater_bowl_parallel
               TARGET test_shallowwater_bowl
               LABELS shallowwater
               TIMEOUT 1500
               CMAKE_GUARD MPI_FOUND
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS   --script fuzzy
                          --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_shallowwater_bowl-reference.vtu
                                  ${CMAKE_CURRENT_BINARY_DIR}/s0002-bowl-parallel-00013.pvtu
                          --zeroThreshold {"velocityY":1e-14,"process rank":100}
                          --command "${MPIEXEC} -np 2 ${CMAKE_CURRENT_BINARY_DIR}/test_shallowwater_bowl -Problem.Name bowl-parallel")
