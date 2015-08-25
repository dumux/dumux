add_dumux_test(test_fluidsystems test_fluidsystems test_fluidsystems.cc
  ${CMAKE_CURRENT_BINARY_DIR}/test_fluidsystems)

#install sources
install(FILES
checkfluidsystem.hh
test_fluidsystems.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/material/fluidsystems)
