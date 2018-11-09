dune_add_test(NAME test_tpfafvgeometry
              SOURCES main.cc
              COMPILE_DEFINITIONS ENABLE_CACHING=false)

dune_add_test(NAME test_tpfafvgeometry_caching
              SOURCES main.cc
              COMPILE_DEFINITIONS ENABLE_CACHING=true)
#install sources
install(FILES
main.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/discretization/cellcentered/tpfa)