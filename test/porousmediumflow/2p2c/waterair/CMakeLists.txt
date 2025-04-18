# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links(FILES params.input)

# non-isothermal tests
dumux_add_test(NAME test_2p2cni_waterair_box
              LABELS porousmediumflow 2p2c  2p2cni
              SOURCES main.cc
              TIMEOUT 3600
              COMPILE_DEFINITIONS TYPETAG=WaterAirBox
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p2cni_waterair_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p2cni_waterair_box-00012.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p2cni_waterair_box params.input -Problem.Name test_2p2cni_waterair_box")

dumux_add_test(NAME test_2p2cni_waterair_buoyancy_box
              LABELS porousmediumflow 2p2c  2p2cni
              TARGET test_2p2cni_waterair_box
              TIMEOUT 3600
              COMPILE_DEFINITIONS TYPETAG=WaterAirBox
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p2cni_waterair_buoyancy_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p2cni_waterair_buoyancy_box-00013.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p2cni_waterair_box params.input -Problem.Name test_2p2cni_waterair_buoyancy_box -TimeLoop.DtInitial 100")

dumux_add_test(NAME test_2p2cni_waterair_tpfa
              LABELS porousmediumflow 2p2c  2p2cni
              SOURCES main.cc
              TIMEOUT 3600
              COMPILE_DEFINITIONS TYPETAG=WaterAirCCTpfa
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p2cni_waterair_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p2cni_waterair_tpfa-00012.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p2cni_waterair_tpfa params.input -Problem.Name test_2p2cni_waterair_tpfa")
