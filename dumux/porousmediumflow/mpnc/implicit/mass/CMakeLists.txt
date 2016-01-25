
#install headers
install(FILES
indicesmass.hh
indicesmasskinetic.hh
localresidualmass.hh
localresidualmasskinetic.hh
volumevariablesmass.hh
volumevariablesmasskinetic.hh
vtkwritermass.hh
vtkwritermasskinetic.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porousmediumflow/mpnc/implicit/mass)
