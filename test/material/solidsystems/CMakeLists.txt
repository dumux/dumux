dune_add_test(SOURCES test_solidsystems.cc)

#install sources
install(FILES
checksolidsystem.hh
test_solidsystems.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/material/solidsystems)
