# build the test for the property system
dune_add_test(SOURCES test_propertysystem.cc)
# build the test for the new property system
dune_add_test(SOURCES test_newpropertysystem.cc)

#install sources
install(FILES
test_propertysystem.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/common/propertysystem)
