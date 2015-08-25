
#install headers
install(FILES
mpncindicesmass.hh
mpncindicesmasskinetic.hh
mpnclocalresidualmass.hh
mpnclocalresidualmasskinetic.hh
mpncvolumevariablesmass.hh
mpncvolumevariablesmasskinetic.hh
mpncvtkwritermass.hh
mpncvtkwritermasskinetic.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/implicit/mpnc/mass)
