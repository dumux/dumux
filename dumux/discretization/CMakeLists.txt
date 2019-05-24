add_subdirectory(box)
add_subdirectory(cellcentered)
add_subdirectory(fem)
add_subdirectory(staggered)
add_subdirectory(projection)

install(FILES
basefvgridgeometry.hh
box.hh
ccmpfa.hh
cctpfa.hh
checkoverlapsize.hh
elementsolution.hh
evalgradients.hh
evalsolution.hh
fluxstencil.hh
functionspacebasis.hh
fvgridvariables.hh
fvproperties.hh
localview.hh
method.hh
scvandscvfiterators.hh
staggered.hh
subcontrolvolumebase.hh
subcontrolvolumefacebase.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/discretization)
