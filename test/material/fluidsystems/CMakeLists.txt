dune_add_test(SOURCES test_fluidsystems.cc
              LABELS unit)

#install sources
install(FILES
checkfluidsystem.hh
test_fluidsystems.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/material/fluidsystems)
