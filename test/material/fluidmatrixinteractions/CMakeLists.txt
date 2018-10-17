add_subdirectory("2p")

dune_add_test(SOURCES test_effectivediffusivitymillingtonquirk.cc
              LABELS unit
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzyData --delimiter " "
                       --files ${CMAKE_SOURCE_DIR}/test/references/effectivediffusivitymillingtonquirk-reference.dat
                               ${CMAKE_CURRENT_BINARY_DIR}/millingtonquirk_d_eff.dat
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_effectivediffusivitymillingtonquirk")

dune_add_test(SOURCES test_effectivediffusivityconstanttortuosity.cc
              LABELS unit
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzyData --delimiter " "
                       --files ${CMAKE_SOURCE_DIR}/test/references/effectivediffusivityconstanttau-reference.dat
                               ${CMAKE_CURRENT_BINARY_DIR}/constant_d_eff.dat
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_effectivediffusivityconstanttortuosity -SpatialParams.Tortuosity 0.3")

#install sources
install(FILES
test_effectivediffusivityconstanttortuosity.cc
test_effectivediffusivitymillingtonquirk.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/material/fluidmatrixinteractions)
