dune_add_test(SOURCES test_immiscibleflash.cc
              LABELS unit)

#install sources
install(FILES
test_immiscibleflash.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/material/immiscibleflash)
