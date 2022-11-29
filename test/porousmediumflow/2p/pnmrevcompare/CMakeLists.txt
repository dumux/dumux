dune_symlink_to_source_files(FILES "params.input")

# using tpfa
add_executable(test_2p_incompressible_tpfa_compare EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_2p_incompressible_tpfa_compare PUBLIC TYPETAG=TwoPIncompressibleBox)

dumux_add_test(NAME test_2p_incompressible_tpfa_compare
              TARGET test_2p_incompressible_tpfa_compare
              LABELS porousmediumflow 2p
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_incompressible_cc-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_incompressible_tpfa-00007.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_incompressible_tpfa params.input -Problem.Name test_2p_incompressible_tpfa")
