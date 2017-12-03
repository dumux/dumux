dune_add_test(NAME test_staggeredfvgeometry
              SOURCES test_staggeredfvgeometry.cc
              COMPILE_DEFINITIONS ENABLE_CACHING=false)

dune_add_test(NAME test_staggeredfvgeometry_caching
              SOURCES test_staggeredfvgeometry.cc
              COMPILE_DEFINITIONS ENABLE_CACHING=true)

#install sources
install(FILES
test_staggeredfvgeometry.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/discretization/staggered)
