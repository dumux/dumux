add_subdirectory(grid)
add_subdirectory(vtk)
add_subdirectory(xml)

install(FILES
adaptivegridrestart.hh
container.hh
defaultiofields.hh
fieldnames.hh
gnuplotinterface.hh
loadsolution.hh
ploteffectivediffusivitymodel.hh
plotmateriallaw.hh
plotmateriallaw3p.hh
plotthermalconductivitymodel.hh
pointcloudvtkwriter.hh
restart.hh
staggeredvtkoutputmodule.hh
vtkfunction.hh
vtkmultiwriter.hh
vtknestedfunction.hh
vtkoutputmodule.hh
vtksequencewriter.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/io)
