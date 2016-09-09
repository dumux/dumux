add_dumux_test(test_container_io test_container_io test_container_io.cc
  ${CMAKE_CURRENT_BINARY_DIR}/test_container_io)

#install sources
install(FILES
test_container_io.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/io/container)
