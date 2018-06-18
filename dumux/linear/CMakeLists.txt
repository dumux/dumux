#install headers
install(FILES
amgbackend.hh
amgparallelhelpers.hh
amgtraits.hh
linearsolveracceptsmultitypematrix.hh
matrixconverter.hh
scotchbackend.hh
seqsolverbackend.hh
solver.hh
vectorexchange.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/linear)
