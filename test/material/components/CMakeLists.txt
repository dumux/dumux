dune_add_test(NAME plot_air
              SOURCES plotproperties.cc
              COMPILE_DEFINITIONS COMPONENT=Dumux::Air<double>
              COMMAND ./plot_air)

dune_add_test(NAME plot_benzene
              SOURCES plotproperties.cc
              COMPILE_DEFINITIONS COMPONENT=Dumux::Benzene<double>
              COMMAND ./plot_benzene false)

dune_add_test(NAME plot_brine
              SOURCES plotproperties.cc
              COMPILE_DEFINITIONS COMPONENT=Dumux::Brine<double>
              COMMAND ./plot_brine false)

dune_add_test(NAME plot_ch4
              SOURCES plotproperties.cc
              COMPILE_DEFINITIONS COMPONENT=Dumux::CH4<double>
              COMMAND ./plot_ch4 false)

dune_add_test(NAME plot_dnapl
              SOURCES plotproperties.cc
              COMPILE_DEFINITIONS COMPONENT=Dumux::DNAPL<double>
              COMMAND ./plot_dnapl false)

dune_add_test(NAME plot_h2
              SOURCES plotproperties.cc
              COMPILE_DEFINITIONS COMPONENT=Dumux::H2<double>
              COMMAND ./plot_h2 false)

dune_add_test(NAME plot_h2o
              SOURCES plotproperties.cc
              COMPILE_DEFINITIONS COMPONENT=Dumux::H2O<double>
              COMMAND ./plot_h2o false)

dune_add_test(NAME plot_heavyoil
              SOURCES plotproperties.cc
              COMPILE_DEFINITIONS COMPONENT=Dumux::HeavyOil<double>
              COMMAND ./plot_heavyoil false)

dune_add_test(NAME plot_lnapl
              SOURCES plotproperties.cc
              COMPILE_DEFINITIONS COMPONENT=Dumux::LNAPL<double>
              COMMAND ./plot_lnapl false)

dune_add_test(NAME plot_mesitylene
              SOURCES plotproperties.cc
              COMPILE_DEFINITIONS COMPONENT=Dumux::Mesitylene<double>
              COMMAND ./plot_mesitylene false)

dune_add_test(NAME plot_n2
              SOURCES plotproperties.cc
              COMPILE_DEFINITIONS COMPONENT=Dumux::N2<double>
              COMMAND ./plot_n2 false)

dune_add_test(NAME plot_o2
              SOURCES plotproperties.cc
              COMPILE_DEFINITIONS COMPONENT=Dumux::O2<double>
              COMMAND ./plot_o2 false)

dune_add_test(NAME plot_simpleco2
              SOURCES plotproperties.cc
              COMPILE_DEFINITIONS COMPONENT=Dumux::SimpleCO2<double>
              COMMAND ./plot_simpleco2 false)

dune_add_test(NAME plot_simpleh2o
              SOURCES plotproperties.cc
              COMPILE_DEFINITIONS COMPONENT=Dumux::SimpleH2O<double>
              COMMAND ./plot_simpleh2o false)

dune_add_test(NAME plot_xylene
              SOURCES plotproperties.cc
              COMPILE_DEFINITIONS COMPONENT=Dumux::Xylene<double>
              COMMAND ./plot_xylene false)
