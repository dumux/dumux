dune_add_test(NAME test_boxfvgeometry
              SOURCES test_boxfvgeometry.cc
              COMPILE_DEFINITIONS ENABLE_CACHING=false)

dune_add_test(NAME test_boxfvgeometry_caching
              SOURCES test_boxfvgeometry.cc
              COMPILE_DEFINITIONS ENABLE_CACHING=true)

#install sources
install(FILES
test_boxfvgeometry.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/discretization/box)
