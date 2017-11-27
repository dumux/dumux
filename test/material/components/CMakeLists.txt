add_dumux_test(plot_air plot_air plotproperties.cc
               plot_air false)
target_compile_definitions(plot_air PUBLIC COMPONENT=Dumux::Air<double>)

add_dumux_test(plot_benzene plot_benzene plotproperties.cc
               plot_benzene false)
target_compile_definitions(plot_benzene PUBLIC COMPONENT=Dumux::Benzene<double>)

add_dumux_test(plot_brine plot_brine plotproperties.cc
               plot_brine false)
target_compile_definitions(plot_brine PUBLIC COMPONENT=Dumux::Brine<double>)

add_dumux_test(plot_ch4 plot_ch4 plotproperties.cc
               plot_ch4 false)
target_compile_definitions(plot_ch4 PUBLIC COMPONENT=Dumux::CH4<double>)

# add_dumux_test(plot_co2 plot_co2 plotproperties.cc
#                plot_co2 false)
# target_compile_definitions(plot_co2 COMPONENT=Dumux::CO2<double>)

add_dumux_test(plot_dnapl plot_dnapl plotproperties.cc
               plot_dnapl false)
target_compile_definitions(plot_dnapl PUBLIC COMPONENT=Dumux::DNAPL<double>)

add_dumux_test(plot_h2 plot_h2 plotproperties.cc
               plot_h2 false)
target_compile_definitions(plot_h2 PUBLIC COMPONENT=Dumux::H2<double>)

add_dumux_test(plot_h2o plot_h2o plotproperties.cc
               plot_h2o false)
target_compile_definitions(plot_h2o PUBLIC COMPONENT=Dumux::H2O<double>)

add_dumux_test(plot_heavyoil plot_heavyoil plotproperties.cc
               plot_heavyoil false)
target_compile_definitions(plot_heavyoil PUBLIC COMPONENT=Dumux::HeavyOil<double>)

add_dumux_test(plot_lnapl plot_lnapl plotproperties.cc
               plot_lnapl false)
target_compile_definitions(plot_lnapl PUBLIC COMPONENT=Dumux::LNAPL<double>)

add_dumux_test(plot_mesitylene plot_mesitylene plotproperties.cc
               plot_mesitylene false)
target_compile_definitions(plot_mesitylene PUBLIC COMPONENT=Dumux::Mesitylene<double>)

add_dumux_test(plot_n2 plot_n2 plotproperties.cc
               plot_n2 false)
target_compile_definitions(plot_n2 PUBLIC COMPONENT=Dumux::N2<double>)

add_dumux_test(plot_o2 plot_o2 plotproperties.cc
               plot_o2 false)
target_compile_definitions(plot_o2 PUBLIC COMPONENT=Dumux::O2<double>)

add_dumux_test(plot_simpleco2 plot_simpleco2 plotproperties.cc
               plot_simpleco2 false)
target_compile_definitions(plot_simpleco2 PUBLIC COMPONENT=Dumux::SimpleCO2<double>)

add_dumux_test(plot_simpleh2o plot_simpleh2o plotproperties.cc
               plot_simpleh2o false)
target_compile_definitions(plot_simpleh2o PUBLIC COMPONENT=Dumux::SimpleH2O<double>)

add_dumux_test(plot_xylene plot_xylene plotproperties.cc
               plot_xylene false)
target_compile_definitions(plot_xylene PUBLIC COMPONENT=Dumux::Xylene<double>)
