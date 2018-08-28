dune_symlink_to_source_files(FILES "tracer.input")

dune_add_test(NAME test_tracer_sequential
              SOURCES test_cctracer.cc
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/tracer-sequential-transport-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/tracer-00010.vtu
                               ${CMAKE_SOURCE_DIR}/test/references/tracer-sequential-pressure-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/1p.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_tracer_sequential tracer.input")

#install sources
install(FILES
1ptestproblem.hh
1ptestspatialparams.hh
tracertestproblem.hh
tracertestspatialparams.hh
test_cctracer.cc
tracer.input
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/porousmediumflow/tracer/1ptracer)
