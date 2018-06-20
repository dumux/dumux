add_subdirectory(box)
add_subdirectory(cellcentered)

install(FILES
couplingmanager.hh
couplingmapper.hh
couplingmapperbase.hh
enrichmenthelper.hh
gmshreader.hh
gridcreator.hh
gridvertexadapter.hh
vertexmapper.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/multidomain/facet)
