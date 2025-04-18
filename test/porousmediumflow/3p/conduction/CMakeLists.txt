# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links(FILES params.input)

dumux_add_test(NAME test_3pni_conduction_box
              LABELS porousmediumflow 3p 3pni
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=ThreePNIConductionBox
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS      --script fuzzy
                            --files ${CMAKE_SOURCE_DIR}/test/references/test_3pni_conduction_box-reference.vtu
                            ${CMAKE_CURRENT_BINARY_DIR}/test_3pni_conduction_box-00005.vtu
                            --command "${CMAKE_CURRENT_BINARY_DIR}/test_3pni_conduction_box params.input -Problem.Name test_3pni_conduction_box"
                            --zeroThreshold {"velocity_aq \(m/s\)_1":1e-8})

dumux_add_test(NAME test_3pni_conduction_tpfa
              LABELS porousmediumflow 3p 3pni
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=ThreePNIConductionCCTpfa
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS      --script fuzzy
                            --files ${CMAKE_SOURCE_DIR}/test/references/test_3pni_conduction_tpfa-reference.vtu
                            ${CMAKE_CURRENT_BINARY_DIR}/test_3pni_conduction_tpfa-00005.vtu
                            --command "${CMAKE_CURRENT_BINARY_DIR}/test_3pni_conduction_tpfa params.input -Problem.Name test_3pni_conduction_tpfa"
                            --zeroThreshold {"velocity_aq \(m/s\)_1":1e-8})
