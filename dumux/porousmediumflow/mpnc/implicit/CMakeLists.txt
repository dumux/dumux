add_subdirectory("diffusion")
add_subdirectory("energy")
add_subdirectory("mass")

#install headers
install(FILES
fluxvariables.hh
indices.hh
localresidual.hh
model.hh
modelkinetic.hh
newtoncontroller.hh
properties.hh
propertieskinetic.hh
propertydefaults.hh
propertydefaultskinetic.hh
volumevariables.hh
volumevariablesia.hh
volumevariablesiakinetic.hh
vtkwritercommon.hh
vtkwriter.hh
vtkwritermodule.hh
velomodelnewtoncontroller.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porousmediumflow/mpnc/implicit)
