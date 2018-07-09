add_subdirectory(box)
add_subdirectory(cellcentered)

install(FILES
codimonegridadapter.hh
couplingmanager.hh
couplingmapper.hh
couplingmapperbase.hh
enrichmenthelper.hh
gmshreader.hh
gridcreator.hh
vertexmapper.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/multidomain/facet)
