# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()

dumux_add_test(NAME test_1pni_conduction_box
              LABELS porousmediumflow 1p 1pni
              SOURCES ../main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePNIConductionBox
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_1pni_conduction_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_1pni_conduction_box-00005.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pni_conduction_box params.input -Problem.Name test_1pni_conduction_box"
                       --zeroThreshold {"velocity_liq \(m/s\)":1e-8})

dumux_add_test(NAME test_1pni_conduction_tpfa
              LABELS porousmediumflow 1p 1pni
              SOURCES ../main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePNIConductionCCTpfa
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_1pni_conduction_cc-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_1pni_conduction_tpfa-00005.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pni_conduction_tpfa params.input -Problem.Name test_1pni_conduction_tpfa"
                       --zeroThreshold {"velocity_liq \(m/s\)":1e-8})

dumux_add_test(NAME test_1pni_conduction_mpfa
              LABELS porousmediumflow 1p 1pni
              SOURCES ../main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePNIConductionCCMpfa
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_1pni_conduction_mpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_1pni_conduction_mpfa-00005.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pni_conduction_mpfa params.input -Problem.Name test_1pni_conduction_mpfa")

# the shared main.cc in the parent folder includes the properties.hh of this folder
foreach(target test_1pni_conduction_box test_1pni_conduction_tpfa test_1pni_conduction_mpfa)
  target_include_directories(${target} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
endforeach()
