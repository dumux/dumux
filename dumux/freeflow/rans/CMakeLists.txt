add_subdirectory("zeroeq")

#install headers
install(FILES
# fluxvariables.hh
# fluxvariablescache.hh
# indices.hh
# localresidual.hh
model.hh
problem.hh
volumevariables.hh
vtkoutputfields.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/freeflow/rans)
