dune_symlink_to_source_files(FILES "test_mcwhorter.input")

dune_add_test(NAME test_mcwhorter
              SOURCES test_mcwhorter_fv.cc
              COMPILE_DEFINITIONS TYPETAG=TwoPTransport
              CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_mcwhorter test_mcwhorter.input")

#install sources
install(FILES
problem.hh
spatialparams.hh
test_transport_fv.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/porousmediumflow/2p/sequential/impes/mcwhorter)
