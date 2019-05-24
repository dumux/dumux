add_subdirectory(boundary)
add_subdirectory(embedded)
add_subdirectory(facet)
add_subdirectory(io)

install(FILES
couplingjacobianpattern.hh
couplingmanager.hh
fvassembler.hh
fvgridgeometry.hh
fvgridvariables.hh
fvproblem.hh
glue.hh
newtonsolver.hh
staggeredcouplingmanager.hh
staggeredtraits.hh
subdomainboxlocalassembler.hh
subdomaincclocalassembler.hh
subdomainstaggeredlocalassembler.hh
traits.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/multidomain)
