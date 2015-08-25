add_dumux_test(test_immiscibleflash test_immiscibleflash test_immiscibleflash.cc
  ${CMAKE_CURRENT_BINARY_DIR}/test_immiscibleflash)

#install sources
install(FILES
test_immiscibleflash.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/material/immiscibleflash)
