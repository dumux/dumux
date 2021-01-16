dune_symlink_to_source_files(FILES "params.input")
dune_symlink_to_source_files(FILES grids)

dumux_add_test(NAME test_2p_cornerpoint
              LABELS porousmediumflow 2p
              SOURCES main.cc
              CMAKE_GUARD HAVE_OPM_GRID
              COMPILE_DEFINITIONS HAVE_ECL_INPUT=1
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p_cornerpoint-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p_cornerpoint-00005.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_cornerpoint params.input -Problem.Name test_2p_cornerpoint")
