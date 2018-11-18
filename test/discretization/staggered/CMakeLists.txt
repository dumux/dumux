dune_add_test(NAME test_staggeredfvgeometry
              SOURCES test_staggeredfvgeometry.cc
              LABELS unit discretization)

dune_add_test(NAME test_staggered_free_flow_geometry
              SOURCES test_staggered_free_flow_geometry.cc
              LABELS unit discretization)

#install sources
install(FILES
test_staggeredfvgeometry.cc
test_staggered_free_flow_geometry.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/discretization/staggered)
