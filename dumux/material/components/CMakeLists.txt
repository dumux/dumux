add_subdirectory("iapws")

#install headers
install(FILES
air.hh
benzene.hh
brine.hh
ch4.hh
co2.hh
co2tables.inc
co2tablereader.hh
component.hh
dnapl.hh
h2.hh
h2o.hh
lnapl.hh
mesitylene.hh
n2.hh
nacl.hh
o2.hh
simpleco2.hh
simpleh2o.hh
tabulatedcomponent.hh
unit.hh
xylene.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/material/components)
