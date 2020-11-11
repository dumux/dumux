add_input_file_links()

add_executable(test_ff_backwardfacingstep EXCLUDE_FROM_ALL main.cc)

dumux_add_test(NAME test_ff_stokes_backwardfacingstep_outflow
              TARGET test_ff_backwardfacingstep
              LABELS freeflow navierstokes
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_channel-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_backwardfacingstep_outflow-00002.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_backwardfacingstep params.input -Problem.OutletCondition Outflow
                             -Problem.Name test_ff_stokes_backwardfacingstep_outflow")
