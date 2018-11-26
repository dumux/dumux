add_subdirectory(boundary)
add_subdirectory(embedded)
add_subdirectory(facet)

install(FILES
couplingjacobianpattern.hh
couplingmanager.hh
fvassembler.hh
newtonsolver.hh
staggeredcouplingmanager.hh
staggeredtraits.hh
subdomainboxlocalassembler.hh
subdomaincclocalassembler.hh
subdomainstaggeredlocalassembler.hh
traits.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/multidomain)
