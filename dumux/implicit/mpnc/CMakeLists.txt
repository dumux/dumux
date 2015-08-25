add_subdirectory("diffusion")
add_subdirectory("energy")
add_subdirectory("mass")

#install headers
install(FILES
mpncfluxvariables.hh
mpncindices.hh
mpnclocalresidual.hh
mpncmodel.hh
mpncmodelkinetic.hh
mpncnewtoncontroller.hh
mpncproperties.hh
mpncpropertieskinetic.hh
mpncpropertydefaults.hh
mpncpropertydefaultskinetic.hh
mpncvolumevariables.hh
mpncvolumevariablesia.hh
mpncvolumevariablesiakinetic.hh
mpncvtkwritercommon.hh
mpncvtkwriter.hh
mpncvtkwritermodule.hh
velomodelnewtoncontroller.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/implicit/mpnc)
