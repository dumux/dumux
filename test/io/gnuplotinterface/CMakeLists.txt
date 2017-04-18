file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
add_dumux_test(test_gnuplotinterface test_gnuplotinterface test_gnuplotinterface.cc
  ${CMAKE_CURRENT_BINARY_DIR}/test_gnuplotinterface)

#install sources
install(FILES
test_gnuplotinterface.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/io/gnuplotinterface)
