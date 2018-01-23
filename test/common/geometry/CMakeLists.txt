dune_add_test(SOURCES test_1d3d_intersection.cc)
dune_add_test(SOURCES test_graham_convex_hull.cc)

#install sources
install(FILES
test_1d3d_intersection.cc
test_graham_convex_hull.cc
writetriangulation.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/common/geometry)

set(CMAKE_BUILD_TYPE RelWithDebInfo)
