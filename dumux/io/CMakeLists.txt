add_subdirectory(grid)

install(FILES
adaptivegridrestart.hh
container.hh
defaultvtkoutputfields.hh
gnuplotinterface.hh
gridcreator.hh
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
