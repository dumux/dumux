# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()

dumux_add_test(NAME test_mpnc_2p2c_comparison_box
              LABELS porousmediumflow mpnc
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=MPNCComparisonBox
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_mpnc_2p2c_comparison_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_mpnc_2p2c_comparison_box-00009.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_mpnc_2p2c_comparison_box params.input -Problem.Name test_mpnc_2p2c_comparison_box")

dumux_add_test(NAME test_mpnc_2p2c_comparison_tpfa
              LABELS porousmediumflow mpnc
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=MPNCComparisonCC
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_mpnc_2p2c_comparison_tpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_mpnc_2p2c_comparison_tpfa-00009.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_mpnc_2p2c_comparison_tpfa params.input -Problem.Name test_mpnc_2p2c_comparison_tpfa")
