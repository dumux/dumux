#install headers
install(FILES
couplingmanager1d3d.hh
couplingmanager2d3d.hh
couplingmanagerbase.hh
integrationpointsource.hh
mixeddimensionglue.hh
pointsourcedata.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/multidomain/embedded)
