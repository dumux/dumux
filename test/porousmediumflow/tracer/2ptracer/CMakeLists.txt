dune_symlink_to_source_files(FILES "tracer_2p.input")

dune_add_test(NAME test_2p_tracer_sequential
              SOURCES test_2p_cctracer.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/tracer-sequential-transport-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/tracer-00010.vtu
                               ${CMAKE_SOURCE_DIR}/test/references/tracer-sequential-pressure-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/1p.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_tracer_sequential tracer_2p.input")

#install sources
install(FILES
2ptestproblem.hh
2ptestspatialparams.hh
2ptracertestproblem.hh
2ptracertestspatialparams.hh
test_2p_cctracer.cc
tracer_2p.input
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/porousmediumflow/tracer/2ptracer)
