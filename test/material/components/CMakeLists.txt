dumux_add_test(SOURCES test_componenttraits.cc
              COMPILE_ONLY
              LABELS unit material)

add_executable(plot_component plotproperties.cc)

dumux_add_test(NAME plot_air
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "Air"
              LABELS unit material)

dumux_add_test(NAME plot_ammonia
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "Ammonia"
              LABELS unit material)

dumux_add_test(NAME plot_benzene
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "Benzene"
              LABELS unit material)

dumux_add_test(NAME plot_brine
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "Brine"
              LABELS unit material)

dumux_add_test(NAME plot_calcite
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "Calcite"
              LABELS unit material)

dumux_add_test(NAME plot_calciumion
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "CalciumIon"
              LABELS unit material)

dumux_add_test(NAME plot_cao
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "CaO"
              LABELS unit material)

dumux_add_test(NAME plot_cao2h2
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "CaO2H2"
              LABELS unit material)

dumux_add_test(NAME plot_carbonateion
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "CarbonateIon"
              LABELS unit material)

dumux_add_test(NAME plot_ch4
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "CH4"
              LABELS unit material)

dumux_add_test(NAME plot_chlorideion
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "ChlorideIon"
              LABELS unit material)

dumux_add_test(NAME plot_glucose
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "Glucose"
              LABELS unit material)

dumux_add_test(NAME plot_granite
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "Granite"
              LABELS unit material)

dumux_add_test(NAME plot_h2
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "H2"
              LABELS unit material)

dumux_add_test(NAME plot_h2o
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "H2O"
              LABELS unit material)

dumux_add_test(NAME plot_heavyoil
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "HeavyOil"
              LABELS unit material)

dumux_add_test(NAME plot_mesitylene
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "Mesitylene"
              LABELS unit material)

dumux_add_test(NAME plot_n2
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "N2"
              LABELS unit material)

dumux_add_test(NAME plot_nacl
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "NaCl"
              LABELS unit material)

dumux_add_test(NAME plot_o2
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "O2"
              LABELS unit material)

dumux_add_test(NAME plot_simpleh2o
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "SimpleH2O"
              LABELS unit material)

dumux_add_test(NAME plot_sodiumion
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "SodiumIon"
              LABELS unit material)

dumux_add_test(NAME plot_trichloroethene
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "Trichloroethene"
              LABELS unit material)

dumux_add_test(NAME plot_urea
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "Urea"
              LABELS unit material)

dumux_add_test(NAME plot_xylene
              TARGET plot_component
              COMMAND ./plot_component
              CMD_ARGS "Xylene"
              LABELS unit material)
