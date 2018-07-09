dune_add_test(NAME test_boundary_equaldim_1p_2p
              SOURCES main.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMAKE_GUARD dune-subgrid_FOUND
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/md-boundary-1p2p-injection-2p-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/2p-00017.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/md-boundary-1p2p-injection-1p-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/1p-00017.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_boundary_equaldim_1p_2p test_boundary_equaldim_1p_2p.input")

dune_symlink_to_source_files(FILES "test_boundary_equaldim_1p_2p.input")
