add_subdirectory("buckleyleverett")
add_subdirectory("mcwhorter")

dune_symlink_to_source_files(FILES "test_impes.input")

dune_add_test(NAME test_impes_fv
              SOURCES test_impes_fv.cc
              COMPILE_DEFINITIONS TYPETAG=TwoPTransport
              CMD_ARGS --command "${CMAKE_CURRENT_BINARY_DIR}/test_impes_fv test_impes.input")

#install sources
install(FILES
problem.hh
spatialparams.hh
test_transport_fv.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/porousmediumflow/2p/sequential/impes)
