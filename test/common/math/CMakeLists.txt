# build the test for the property system
dune_add_test(SOURCES test_vtmv.cc)

#install sources
install(FILES
test_vtmv.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/common/math)