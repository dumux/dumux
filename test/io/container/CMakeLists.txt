dune_add_test(SOURCES test_container_io.cc
              LABELS unit)

#install sources
install(FILES
test_container_io.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/io/container)
