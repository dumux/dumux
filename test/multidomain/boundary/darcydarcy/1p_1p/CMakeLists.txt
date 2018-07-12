dune_add_test(NAME test_boundary_equaldim_1p_1p_half
              SOURCES main.cc
              COMPILE_DEFINITIONS DOMAINSPLIT=0
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/1ptestcc-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/md_1p_half_combined.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_boundary_equaldim_1p_1p_half test_boundary_equaldim_1p_1p.input \
                                   -Problem.Name md_1p_half")

dune_add_test(NAME test_boundary_equaldim_1p_1p_lens
              SOURCES main.cc
              COMPILE_DEFINITIONS DOMAINSPLIT=1
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMAKE_GUARD dune-subgrid_FOUND
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/1ptestcc-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/md_1p_lens_combined.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_boundary_equaldim_1p_1p_lens test_boundary_equaldim_1p_1p.input \
                                   -Problem.Name md_1p_lens")

dune_symlink_to_source_files(FILES "test_boundary_equaldim_1p_1p.input")
