add_subdirectory(grid)
add_subdirectory(vtk)
add_subdirectory(xml)

install(FILES
adaptivegridrestart.hh
container.hh
defaultiofields.hh
gnuplotinterface.hh
loadsolution.hh
name.hh
ploteffectivediffusivitymodel.hh
plotmateriallaw.hh
plotmateriallaw3p.hh
plotthermalconductivitymodel.hh
pointcloudvtkwriter.hh
restart.hh
staggeredvtkoutputmodule.hh
velocityoutput.hh
vtkfunction.hh
vtkmultiwriter.hh
vtknestedfunction.hh
vtkoutputmodule.hh
vtksequencewriter.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/io)
