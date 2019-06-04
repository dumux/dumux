dune_symlink_to_source_files(FILES "params.input")

# using box
dune_add_test(NAME test_2p_incompressible_box
              LABELS porousmediumflow 2p
              SOURCES main.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_incompressible_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_incompressible_box-00007.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_incompressible_box params.input -Problem.Name test_2p_incompressible_box -Grid.Cells \"48 32\"")
