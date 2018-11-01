dune_add_test(NAME test_staggeredfvgeometry
              SOURCES main.cc)

#install sources
install(FILES
main.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/discretization/staggered/fv)
