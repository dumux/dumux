add_input_file_links()

# isothermal tests
dune_add_test(SOURCES test_cc2pff.cc)
dune_add_test(SOURCES test_constvel.cc
              COMPILE_DEFINITIONS TYPETAG=ConstVelCCProblem)
dune_add_test(SOURCES test2_constvel.cc)

#install sources
install(FILES
lensproblem.hh
lensspatialparams.hh
test_cc2pff.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/2pff)
