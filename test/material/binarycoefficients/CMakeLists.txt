dune_add_test(SOURCES test_fullermethod.cc)

#install sources
install(FILES
test_fullermethod.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/material/binarycoefficient)
