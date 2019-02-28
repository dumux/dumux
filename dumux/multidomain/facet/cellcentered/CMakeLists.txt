add_subdirectory(mpfa)
add_subdirectory(tpfa)

install(FILES
localresidual.hh
upwindscheme.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/multidomain/facet/cellcentered)
