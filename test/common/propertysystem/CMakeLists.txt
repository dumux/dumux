# build the test for the property system
dune_add_test(SOURCES test_propertysystem.cc
              LABELS unit)
# build the test for the property system using the old macros (deprecated)
dune_add_test(SOURCES test_propertysystem_macros.cc
              LABELS unit)
