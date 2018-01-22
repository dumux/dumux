dune_add_test(SOURCES test_1d3d_intersection.cc)

#install sources
install(FILES
test_1d3d_intersection.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/common/geometry)

set(CMAKE_BUILD_TYPE Release)
