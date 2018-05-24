add_subdirectory(box)
add_subdirectory(cellcentered)

#install headers
install(FILES
gmshreader.hh
gridcreator.hh
gridvertexadapter.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/multidomain/facet)
