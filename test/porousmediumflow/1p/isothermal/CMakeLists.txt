# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()

# isothermal tests
dumux_add_test(NAME test_1p_tpfa
              LABELS porousmediumflow 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePTestCCTpfa
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_cc-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p_tpfa-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_tpfa params.input -Problem.Name test_1p_tpfa")

dumux_add_test(NAME test_1p_mpfa
              LABELS porousmediumflow 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePTestCCMpfa
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_cc-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p_mpfa-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_mpfa params.input -Problem.Name test_1p_mpfa")

dumux_add_test(NAME test_1p_box
              LABELS porousmediumflow 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePTestBox
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p_box-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_box params.input -Problem.Name test_1p_box")

dumux_add_test(NAME test_1p_forchheimer_tpfa
              LABELS porousmediumflow 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePTestCCTpfa FORCHHEIMER=1
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_forchheimer_tpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p_forchheimer_tpfa-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_forchheimer_tpfa params_forchheimer.input -Problem.Name test_1p_forchheimer_tpfa"
                        --zeroThreshold {"velocity_liq \(m/s\)":1e-12})

dumux_add_test(NAME test_1p_forchheimer_box
              LABELS porousmediumflow 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePTestBox FORCHHEIMER=1
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_forchheimer_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_1p_forchheimer_box-00001.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_forchheimer_box params_forchheimer.input -Problem.Name test_1p_forchheimer_box"
                        --zeroThreshold {"velocity_liq \(m/s\)":1e-12})

# a gstat test (because it's a random permeability field we can't test against a reference solution)
dune_symlink_to_source_files(FILES control.gstat)
dumux_add_test(NAME test_1p_gstat
              LABELS porousmediumflow 1p
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePTestCCTpfa
              CMAKE_GUARD DUMUX_HAVE_GSTAT
              COMMAND ./test_1p_gstat
              CMD_ARGS params_gstat.input)
