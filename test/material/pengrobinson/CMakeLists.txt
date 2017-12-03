dune_add_test(SOURCES test_pengrobinson.cc)

#install sources
install(FILES
test_pengrobinson.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/material/pengrobinson)
