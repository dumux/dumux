dune_add_test(SOURCES test_boxfvgeometry.cc)

#install sources
install(FILES
test_boxfvgeometry.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/discretization/box)
