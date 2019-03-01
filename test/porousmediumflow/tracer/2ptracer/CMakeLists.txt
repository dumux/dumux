dune_symlink_to_source_files(FILES "params.input")

dune_add_test(NAME test_2p_tracer_sequential
              SOURCES main.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/tracer-sequential-transport-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/tracer-00010.vtu
                               ${CMAKE_SOURCE_DIR}/test/references/tracer-sequential-pressure-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/1p.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_tracer_sequential params.input")

#install sources
install(FILES
problem_2p.hh
problem_tracer.hh
spatialparams_2p.hh
spatialparams_tracer.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/porousmediumflow/tracer/2ptracer)
