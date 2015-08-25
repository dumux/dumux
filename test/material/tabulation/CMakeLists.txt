add_dumux_test(test_tabulation test_tabulation test_tabulation.cc
  ${CMAKE_CURRENT_BINARY_DIR}/test_tabulation)

#install sources
install(FILES
test_tabulation.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/material/tabulation)
