dune_symlink_to_source_files(FILES "test_2pnc.input")

dune_add_test(NAME test_2pnc
              SOURCES test_2pnc.cc
              COMPILE_DEFINITIONS TYPETAG=TwoPNCTestProblemTpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/tracer-sequential-transport-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/tracer-00010.vtu
                               ${CMAKE_SOURCE_DIR}/test/references/tracer-sequential-pressure-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/1p.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pnc tracer_2p.input")

#install sources
install(FILES
2ptestproblem.hh
2ptestspatialparams.hh
2ptestspatialparams_randomfield.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/porousmediumflow/tracer/2ptracer/2pnccomparison)
