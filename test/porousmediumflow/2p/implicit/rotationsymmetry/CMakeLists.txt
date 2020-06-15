dune_symlink_to_source_files(FILES "params.input" "dome.msh")

dumux_add_test(NAME test_2p_rotationsymmetry_dome
               SOURCES main.cc
               LABELS porousmediumflow 2p
               CMAKE_GUARD dune-alugrid_FOUND
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_rotationsymmetry_dome-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_2p_rotationsymmetry_dome-00002.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_rotationsymmetry_dome -TimeLoop.TEnd 1e6")
