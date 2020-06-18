dune_symlink_to_source_files(FILES "params.input")

dune_add_test(NAME example_1p_rotationsymmetry
              LABELS porousmediumflow 1p example
              SOURCES main.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/example_rotationsymmetry-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/example_rotationsymmetry-00000.vtp
                        --command "${CMAKE_CURRENT_BINARY_DIR}/example_1p_rotationsymmetry params.input -Problem.Name example_rotationsymmetry")
