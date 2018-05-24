add_subdirectory(tpfa)

#install headers
install(FILES
localresidual.hh
upwindscheme.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/multidomain/facet/cellcentered)
