add_input_file_links()

# isothermal tests
add_dumux_test(test_cc2pff test_cc2pff test_cc2pff.cc
               ${CMAKE_CURRENT_BINARY_DIR}/test_cc2pff)

#install sources
install(FILES
lensproblem.hh
lensspatialparams.hh
test_cc2pff.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/2pff)
