# build the test for the property system
dune_add_test(SOURCES test_math.cc
              LABELS unit)

#install sources
install(FILES
test_math.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/common/math)
