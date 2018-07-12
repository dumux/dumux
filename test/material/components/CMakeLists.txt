dune_add_test(SOURCES test_componenttraits.cc
              COMPILE_ONLY)

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

dune_add_test(NAME plot_calcite
            TARGET plot_component
            COMMAND ./plot_component
            CMD_ARGS "Calcite")

dune_add_test(NAME plot_calciumion
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "CalciumIon")

dune_add_test(NAME plot_cao
            TARGET plot_component
            COMMAND ./plot_component
            CMD_ARGS "CaO")

dune_add_test(NAME plot_cao2h2
            TARGET plot_component
            COMMAND ./plot_component
            CMD_ARGS "CaO2H2")

dune_add_test(NAME plot_carbonateion
            TARGET plot_component
            COMMAND ./plot_component
            CMD_ARGS "CarbonateIon")

dune_add_test(NAME plot_ch4
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "CH4")

dune_add_test(NAME plot_granite
            TARGET plot_component
            COMMAND ./plot_component
            CMD_ARGS "Granite")

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

dune_add_test(NAME plot_nacl
            TARGET plot_component
            COMMAND ./plot_component
            CMD_ARGS "NaCl")

dune_add_test(NAME plot_o2
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "O2")

dune_add_test(NAME plot_simpleh2o
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "SimpleH2O")

dune_add_test(NAME plot_trichloroethene
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "Trichloroethene")

dune_add_test(NAME plot_xylene
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "Xylene")
