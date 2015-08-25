add_dumux_test(test_ncpflash test_ncpflash test_ncpflash.cc
  ${CMAKE_CURRENT_BINARY_DIR}/test_ncpflash)

#install sources
install(FILES
test_ncpflash.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/material/ncpflash)
