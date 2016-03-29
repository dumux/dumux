
#install headers
install(FILES
adaptiveproperties.hh
problem.hh
properties.hh
celldataadaptive.hh
celldata.hh
celldatamultiphysics.hh
fluxdata.hh
fv2dpressureadaptive.hh
fv2dtransportadaptive.hh
fv3dpressureadaptive.hh
fv3dtransportadaptive.hh
fvmpfal3dinteractionvolumecontaineradaptive.hh
fvpressure.hh
fvpressuremultiphysics.hh
fvpressurecompositional.hh
fvtransport.hh
fvtransportmultiphysics.hh
variableclassadaptive.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porousmediumflow/2p2c/sequential)
