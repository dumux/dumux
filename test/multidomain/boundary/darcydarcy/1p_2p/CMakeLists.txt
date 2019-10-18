dumux_add_test(NAME test_md_boundary_darcy1p_darcy2p
              LABELS multidomain multidomain_boundary darcydarcy 1p 2p
              SOURCES main.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMAKE_GUARD dune-subgrid_FOUND
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p_darcy2p_2p-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p_darcy2p_2p-00017.vtu
                                ${CMAKE_SOURCE_DIR}/test/references/test_md_boundary_darcy1p_darcy2p_1p-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p_darcy2p_1p-00017.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p_darcy2p params.input
                        -Vtk.OutputName test_md_boundary_darcy1p_darcy2p")

dune_symlink_to_source_files(FILES "params.input")
