dumux_add_test(NAME test_md_boundary_darcy1p_darcy1p_half
              LABELS multidomain multidomain_boundary darcydarcy 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS DOMAINSPLIT=0
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_cc-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p_darcy1p_half_combined.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p_darcy1p_half params.input \
                                   -Vtk.OutputName test_md_boundary_darcy1p_darcy1p_half")

dumux_add_test(NAME test_md_boundary_darcy1p_darcy1p_lens
              LABELS multidomain multidomain_boundary darcydarcy 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS DOMAINSPLIT=1
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMAKE_GUARD dune-subgrid_FOUND
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_cc-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p_darcy1p_lens_combined.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_boundary_darcy1p_darcy1p_lens params.input \
                                    -Vtk.OutputName test_md_boundary_darcy1p_darcy1p_lens")

dune_symlink_to_source_files(FILES "params.input")
