dune_add_test(SOURCES test_spline.cc
              LABELS unit)

#install sources
install(FILES
test_spline.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/common/spline)
