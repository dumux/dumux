add_subdirectory(box)
add_subdirectory(cellcentered)
add_subdirectory(staggered)

install(FILES
basefvgridgeometry.hh
checkoverlapsize.hh
elementsolution.hh
evalgradients.hh
evalsolution.hh
fluxstencil.hh
fvgridvariables.hh
fvproperties.hh
localview.hh
method.hh
scvandscvfiterators.hh
subcontrolvolumebase.hh
subcontrolvolumefacebase.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/discretization)
