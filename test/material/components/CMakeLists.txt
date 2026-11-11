# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(SOURCES test_componenttraits.cc
              COMPILE_ONLY
              LABELS unit material)

dumux_add_test_target(NAME test_plot_component
                      SOURCES test_plot_component.cc)

dumux_add_test(NAME test_h2o_simpleh2o_consistency
               SOURCES test_h2o_simpleh2o_consistency.cc
               LABELS unit material)

dumux_add_test(NAME test_plot_air
              TARGET test_plot_component
              COMMAND ./test_plot_component
              CMD_ARGS "Air"
              LABELS unit material)

dumux_add_test(NAME test_plot_ammonia
              TARGET test_plot_component
              COMMAND ./test_plot_component
              CMD_ARGS "Ammonia"
              LABELS unit material)

dumux_add_test(NAME test_plot_benzene
              TARGET test_plot_component
              COMMAND ./test_plot_component
              CMD_ARGS "Benzene"
              LABELS unit material)

dumux_add_test(NAME test_plot_brine
              TARGET test_plot_component
              COMMAND ./test_plot_component
              CMD_ARGS "Brine"
              LABELS unit material)

dumux_add_test(NAME test_plot_calcite
              TARGET test_plot_component
              COMMAND ./test_plot_component
              CMD_ARGS "Calcite"
              LABELS unit material)

dumux_add_test(NAME test_plot_calciumion
              TARGET test_plot_component
              COMMAND ./test_plot_component
              CMD_ARGS "CalciumIon"
              LABELS unit material)

dumux_add_test(NAME test_plot_cao
              TARGET test_plot_component
              COMMAND ./test_plot_component
              CMD_ARGS "CaO"
              LABELS unit material)

dumux_add_test(NAME test_plot_cao2h2
              TARGET test_plot_component
              COMMAND ./test_plot_component
              CMD_ARGS "CaO2H2"
              LABELS unit material)

dumux_add_test(NAME test_plot_carbonateion
              TARGET test_plot_component
              COMMAND ./test_plot_component
              CMD_ARGS "CarbonateIon"
              LABELS unit material)

dumux_add_test(NAME test_plot_ch4
              TARGET test_plot_component
              COMMAND ./test_plot_component
              CMD_ARGS "CH4"
              LABELS unit material)

dumux_add_test(NAME test_plot_constantcomponent
              TARGET test_plot_component
              COMMAND ./test_plot_component
              CMD_ARGS "Constant"
              LABELS unit material)

dumux_add_test(NAME test_plot_chlorideion
              TARGET test_plot_component
              COMMAND ./test_plot_component
              CMD_ARGS "ChlorideIon"
              LABELS unit material)

dumux_add_test(NAME test_plot_glucose
              TARGET test_plot_component
              COMMAND ./test_plot_component
              CMD_ARGS "Glucose"
              LABELS unit material)

dumux_add_test(NAME test_plot_granite
              TARGET test_plot_component
              COMMAND ./test_plot_component
              CMD_ARGS "Granite"
              LABELS unit material)

dumux_add_test(NAME test_plot_h2
              TARGET test_plot_component
              COMMAND ./test_plot_component
              CMD_ARGS "H2"
              LABELS unit material)

dumux_add_test(NAME test_plot_h2o
              TARGET test_plot_component
              COMMAND ./test_plot_component
              CMD_ARGS "H2O"
              LABELS unit material)

dumux_add_test(NAME test_plot_heavyoil
              TARGET test_plot_component
              COMMAND ./test_plot_component
              CMD_ARGS "HeavyOil"
              LABELS unit material)

dumux_add_test(NAME test_plot_mesitylene
              TARGET test_plot_component
              COMMAND ./test_plot_component
              CMD_ARGS "Mesitylene"
              LABELS unit material)

dumux_add_test(NAME test_plot_n2
              TARGET test_plot_component
              COMMAND ./test_plot_component
              CMD_ARGS "N2"
              LABELS unit material)

dumux_add_test(NAME test_plot_nacl
              TARGET test_plot_component
              COMMAND ./test_plot_component
              CMD_ARGS "NaCl"
              LABELS unit material)

dumux_add_test(NAME test_plot_o2
              TARGET test_plot_component
              COMMAND ./test_plot_component
              CMD_ARGS "O2"
              LABELS unit material)

dumux_add_test(NAME test_plot_simpleco2
              TARGET test_plot_component
              COMMAND ./test_plot_component
              CMD_ARGS "SimpleCO2"
              LABELS unit material)

dumux_add_test(NAME test_plot_simpleh2o
              TARGET test_plot_component
              COMMAND ./test_plot_component
              CMD_ARGS "SimpleH2O"
              LABELS unit material)

dumux_add_test(NAME test_plot_sodiumion
              TARGET test_plot_component
              COMMAND ./test_plot_component
              CMD_ARGS "SodiumIon"
              LABELS unit material)

dumux_add_test(NAME test_plot_trichloroethene
              TARGET test_plot_component
              COMMAND ./test_plot_component
              CMD_ARGS "Trichloroethene"
              LABELS unit material)

dumux_add_test(NAME test_plot_urea
              TARGET test_plot_component
              COMMAND ./test_plot_component
              CMD_ARGS "Urea"
              LABELS unit material)

dumux_add_test(NAME test_plot_xylene
              TARGET test_plot_component
              COMMAND ./test_plot_component
              CMD_ARGS "Xylene"
              LABELS unit material)
