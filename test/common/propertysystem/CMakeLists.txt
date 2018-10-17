# build the test for the property system
dune_add_test(SOURCES test_propertysystem.cc
              LABELS unit)
# build the test for the property system using the old macros (deprecated)
dune_add_test(SOURCES test_propertysystem_macros.cc
              LABELS unit)

#install sources
install(FILES
test_propertysystem.cc
test_propertysystem_macros.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/common/propertysystem)
