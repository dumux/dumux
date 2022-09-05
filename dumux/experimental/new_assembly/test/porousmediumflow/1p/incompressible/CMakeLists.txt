dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME test_new_assembly_1p_incompressible_tpfa
               SOURCES main.cc
               LABELS porousmediumflow 1p new_assembly
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS  --script fuzzy
                         --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_cc-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_1p_incompressible_tpfa-00001.vtu
                         --command "${CMAKE_CURRENT_BINARY_DIR}/test_new_assembly_1p_incompressible_tpfa params.input -Problem.Name test_1p_incompressible_tpfa")
