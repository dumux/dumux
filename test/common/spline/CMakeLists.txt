add_dumux_test(test_spline test_spline test_spline.cc
  ${CMAKE_CURRENT_BINARY_DIR}/test_spline)

#install sources
install(FILES
test_spline.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/common/spline)
