dune_add_test(SOURCES test_tabulation.cc
              LABELS unit)

#install sources
install(FILES
test_tabulation.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/material/tabulation)
