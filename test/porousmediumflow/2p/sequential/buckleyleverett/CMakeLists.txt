dune_symlink_to_source_files(FILES "test_buckleyleverett.input")

dune_add_test(NAME test_buckleyleverett
              SOURCES test_buckleyleverett_fv.cc
              COMPILE_DEFINITIONS TYPETAG=TwoPTransport
              CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_buckleyleverett test_buckleyleverett.input")

#install sources
install(FILES
problem.hh
spatialparams.hh
test_transport_fv.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/porousmediumflow/2p/sequential/buckleyleverett)
