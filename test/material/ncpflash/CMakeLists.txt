dune_add_test(SOURCES test_ncpflash.cc
              LABELS unit)

#install sources
install(FILES
test_ncpflash.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/material/ncpflash)
