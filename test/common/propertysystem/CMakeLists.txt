# build the test for the property system
dumux_add_test(SOURCES test_propertysystem.cc
              LABELS unit)

dumux_add_test(SOURCES test_propertysystem_aliases.cc
              LABELS unit)
