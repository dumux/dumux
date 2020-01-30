add_subdirectory(2p)

dumux_add_test(SOURCES test_effectivediffusivity.cc
              LABELS unit material
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzyData --delimiter " "
                       --files ${CMAKE_SOURCE_DIR}/test/references/effectivediffusivityconstanttau-reference.dat
                               ${CMAKE_CURRENT_BINARY_DIR}/constant_d_eff.dat
                               ${CMAKE_SOURCE_DIR}/test/references/effectivediffusivitymillingtonquirk-reference.dat
                               ${CMAKE_CURRENT_BINARY_DIR}/millingtonquirk_d_eff.dat
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_effectivediffusivity -SpatialParams.Tortuosity 0.3")
