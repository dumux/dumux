dumux_add_test(NAME test_ff_stokes2c_maxwellstefan
              LABELS freeflow navierstokes
              SOURCES main.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes2c_maxwellstefan-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes2c_maxwellstefan-00005.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes2c_maxwellstefan params.input
                             -Problem.Name test_ff_stokes2c_maxwellstefan")

dune_symlink_to_source_files(FILES "params.input")
