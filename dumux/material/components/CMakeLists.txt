add_subdirectory("iapws")

#install headers
install(FILES
air.hh
benzene.hh
brine.hh
CaO.hh
CaOtest.hh
CaO2H2.hh
ch4.hh
co2.hh
co2tablereader.hh
co2tables.inc
component.hh
constant.hh
dnapl.hh
h2.hh
h2o.hh
heavyoil.hh
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
