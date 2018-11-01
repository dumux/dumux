dune_add_test(NAME test_tpfafvgeometry_nonconforming
              SOURCES main.cc
              COMPILE_DEFINITIONS ENABLE_CACHING=false
              CMAKE_GUARD dune-alugrid_FOUND)

dune_add_test(NAME test_cachedtpfafvgeometry_nonconforming
              SOURCES main.cc
              COMPILE_DEFINITIONS ENABLE_CACHING=true
              CMAKE_GUARD dune-alugrid_FOUND)
#install sources
install(FILES
main.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/discretization/cellcentered/tpfanonconforming)