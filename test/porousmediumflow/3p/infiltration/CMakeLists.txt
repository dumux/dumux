# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links(FILES params.input)

# isothermal tests
dumux_add_test(NAME test_3p_infiltration_box
              LABELS porousmediumflow 3p
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=InfiltrationThreePBox
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS      --script fuzzy
                            --files ${CMAKE_SOURCE_DIR}/test/references/test_3p_infiltration_box-reference.vtu
                                    ${CMAKE_CURRENT_BINARY_DIR}/test_3p_infiltration_box-00008.vtu
                            --command "${CMAKE_CURRENT_BINARY_DIR}/test_3p_infiltration_box params.input -Problem.Name test_3p_infiltration_box")

dumux_add_test(NAME test_3p_infiltration_tpfa
              LABELS porousmediumflow 3p
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=InfiltrationThreePCCTpfa
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS      --script fuzzy
                            --files ${CMAKE_SOURCE_DIR}/test/references/test_3p_infiltration_tpfa-reference.vtu
                                    ${CMAKE_CURRENT_BINARY_DIR}/test_3p_infiltration_tpfa-00008.vtu
                            --command "${CMAKE_CURRENT_BINARY_DIR}/test_3p_infiltration_tpfa params.input -Problem.Name test_3p_infiltration_tpfa")
