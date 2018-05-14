add_subdirectory("poromechanics")

#install headers
install(FILES
traits.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/multidomain)
