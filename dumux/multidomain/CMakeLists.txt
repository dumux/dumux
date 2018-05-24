add_subdirectory(poromechanics)
add_subdirectory(embedded)
add_subdirectory(facet)

#install headers
install(FILES
traits.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/multidomain)
