# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()

dumux_add_test(NAME test_2p2c_nonequilibrium_box
              LABELS porousmediumflow 2p2c
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TwoPTwoCChemicalNonequilibriumBox
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p2c_nonequilibrium_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p2c_nonequilibrium_box-00005.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p2c_nonequilibrium_box params.input -Problem.Name test_2p2c_nonequilibrium_box")

dumux_add_test(NAME test_2p2c_nonequilibrium_tpfa
              LABELS porousmediumflow 2p2c
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TwoPTwoCChemicalNonequilibriumCC
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2p2c_nonequilibrium_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p2c_nonequilibrium_tpfa-00005.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p2c_nonequilibrium_tpfa params.input -Problem.Name test_2p2c_nonequilibrium_tpfa")
