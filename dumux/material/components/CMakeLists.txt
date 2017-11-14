add_subdirectory("iapws")

add_dumux_test(plot_air plot_air plotproperties.cc
               plot_air false)
target_compile_definitions(plot_air PUBLIC "AIR=1")

add_dumux_test(plot_benzene plot_benzene plotproperties.cc
               plot_benzene false)
target_compile_definitions(plot_benzene PUBLIC "BENZENE=1")

add_dumux_test(plot_brine plot_brine plotproperties.cc
               plot_brine false)
target_compile_definitions(plot_brine PUBLIC "BRINE=1")

add_dumux_test(plot_ch4 plot_ch4 plotproperties.cc
               plot_ch4 false)
target_compile_definitions(plot_ch4 PUBLIC "METHANE=1")

# add_dumux_test(plot_co2 plot_co2 plotproperties.cc
#                plot_co2 false)
# target_compile_definitions(plot_co2 PUBLIC "CARBONDIOXIDE=1")

add_dumux_test(plot_dnapl plot_dnapl plotproperties.cc
               plot_dnapl false)
target_compile_definitions(plot_dnapl PUBLIC "DNAPL_TCE=1")

add_dumux_test(plot_h2 plot_h2 plotproperties.cc
               plot_h2 false)
target_compile_definitions(plot_h2 PUBLIC "HYDROGEN=1")

add_dumux_test(plot_h2o plot_h2o plotproperties.cc
               plot_h2o false)
target_compile_definitions(plot_h2o PUBLIC "WATER=1")

add_dumux_test(plot_heavyoil plot_heavyoil plotproperties.cc
               plot_heavyoil false)
target_compile_definitions(plot_heavyoil PUBLIC "HEAVYOIL=1")

add_dumux_test(plot_lnapl plot_lnapl plotproperties.cc
               plot_lnapl false)
target_compile_definitions(plot_lnapl PUBLIC "LNAPL_OIL=1")

add_dumux_test(plot_mesitylene plot_mesitylene plotproperties.cc
               plot_mesitylene false)
target_compile_definitions(plot_mesitylene PUBLIC "MESITYLENE=1")

add_dumux_test(plot_n2 plot_n2 plotproperties.cc
               plot_n2 false)
target_compile_definitions(plot_n2 PUBLIC "NITROGEN=1")

add_dumux_test(plot_o2 plot_o2 plotproperties.cc
               plot_o2 false)
target_compile_definitions(plot_o2 PUBLIC "OXYGEN=1")

add_dumux_test(plot_simpleco2 plot_simpleco2 plotproperties.cc
               plot_simpleco2 false)
target_compile_definitions(plot_simpleco2 PUBLIC "SIMPLE_CARBONDIOXIDE=1")

add_dumux_test(plot_simpleh2o plot_simpleh2o plotproperties.cc
               plot_simpleh2o false)
target_compile_definitions(plot_simpleh2o PUBLIC "SIMPLE_WATER=1")

add_dumux_test(plot_xylene plot_xylene plotproperties.cc
               plot_xylene false)
target_compile_definitions(plot_xylene PUBLIC "XYLENE=1")

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
nullcomponent.hh
o2.hh
simpleco2.hh
simpleh2o.hh
tabulatedcomponent.hh
unit.hh
xylene.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/material/components)
