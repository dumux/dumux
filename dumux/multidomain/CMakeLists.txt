add_subdirectory(poromechanics)
add_subdirectory(embedded)

#install headers
install(FILES
traits.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/multidomain)
