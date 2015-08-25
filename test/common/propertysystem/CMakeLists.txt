# build the test for the property system
add_dumux_test(test_propertysystem test_propertysystem test_propertysystem.cc
  ${CMAKE_CURRENT_BINARY_DIR}/test_propertysystem)

#install sources
install(FILES
test_propertysystem.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/common/propertysystem)
