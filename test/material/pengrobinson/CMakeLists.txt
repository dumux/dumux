add_dumux_test(test_pengrobinson test_pengrobinson test_pengrobinson.cc
  ${CMAKE_CURRENT_BINARY_DIR}/test_pengrobinson)

#install sources
install(FILES
test_pengrobinson.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/material/pengrobinson)
