dune_symlink_to_source_files(FILES "test_2p_cornerpoint.input")
dune_symlink_to_source_files(FILES grids)

dune_add_test(NAME test_2p_cornerpoint
              SOURCES test_2p_cornerpoint.cc
              CMAKE_GUARD HAVE_OPM_GRID
              COMPILE_DEFINITIONS HAVE_ECL_INPUT=1
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/cc2pcornerpoint-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/cc2pcornerpoint-00005.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p_cornerpoint -Problem.Name cc2pcornerpoint")

set(CMAKE_BUILD_TYPE Release)

#install sources
install(FILES
problem.hh
spatialparams.hh
test_2p_cornerpoint.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/2p/cornerpoint)
