file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/output")
dune_add_test(SOURCES test_gnuplotinterface.cc
              LABELS unit)

#install sources
install(FILES
test_gnuplotinterface.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/io/gnuplotinterface)
