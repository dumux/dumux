
#install headers
install(FILES
fluxvariablesenergy.hh
fluxvariablesenergykinetic.hh
indicesenergy.hh
indicesenergykinetic.hh
localresidualenergy.hh
localresidualenergykinetic.hh
volumevariablesenergy.hh
volumevariablesenergykinetic.hh
vtkwriterenergy.hh
vtkwriterenergykinetic.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porousmediumflow/mpnc/implicit/energy)
