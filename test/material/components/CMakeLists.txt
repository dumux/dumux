add_executable(plot_component plotproperties.cc)

dune_add_test(NAME plot_air
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "Air")

dune_add_test(NAME plot_benzene
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "Benzene")

dune_add_test(NAME plot_brine
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "Brine")

dune_add_test(NAME plot_ch4
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "CH4")

dune_add_test(NAME plot_dnapl_tce
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "DNAPL_TCE")

dune_add_test(NAME plot_h2
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "H2")

dune_add_test(NAME plot_h2o
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "H2O")

dune_add_test(NAME plot_heavyoil
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "HeavyOil")

dune_add_test(NAME plot_mesitylene
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "Mesitylene")

dune_add_test(NAME plot_n2
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "N2")

dune_add_test(NAME plot_o2
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "O2")

dune_add_test(NAME plot_simpleco2
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "SimpleCO2")

dune_add_test(NAME plot_simpleh2o
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "SimpleH2O")

dune_add_test(NAME plot_xylene
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "Xylene")
