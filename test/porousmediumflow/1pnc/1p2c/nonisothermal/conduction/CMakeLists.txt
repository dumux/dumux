# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

dune_symlink_to_source_files(FILES "params.input")

# non-isothermal tests
# conduction
dumux_add_test(NAME test_1p2cni_conduction_box
              LABELS porousmediumflow 1p2c 1p2cni 1pnc
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePTwoCNIConductionBox
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p2cni_conduction_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p2cni_conduction_box-00005.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p2cni_conduction_box params.input -Problem.Name test_1p2cni_conduction_box"
                        --zeroThreshold {"velocity_liq \(m/s\)":1e-9})

dumux_add_test(NAME test_1p2cni_conduction_tpfa
              LABELS porousmediumflow 1p2c 1p2cni 1pnc
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePTwoCNIConductionCCTpfa
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p2cni_conduction_tpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p2cni_conduction_tpfa-00005.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p2cni_conduction_tpfa params.input -Problem.Name test_1p2cni_conduction_tpfa"
                        --zeroThreshold {"velocity_liq \(m/s\)":1e-9})

dumux_add_test(NAME test_1p2cni_conduction_mpfa
              LABELS porousmediumflow 1p2c 1p2cni 1pnc
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePTwoCNIConductionCCMpfa
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p2cni_conduction_tpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p2cni_conduction_mpfa-00005.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p2cni_conduction_mpfa params.input -Problem.Name test_1p2cni_conduction_mpfa"
                        --zeroThreshold {"velocity_liq \(m/s\)":1e-9})
