# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()

add_executable(test_1pni_convection_box EXCLUDE_FROM_ALL ../main.cc)
target_compile_definitions(test_1pni_convection_box PUBLIC TYPETAG=OnePNIConvectionBox)

dumux_add_test(NAME test_1pni_convection_box
              TARGET test_1pni_convection_box
              LABELS porousmediumflow 1p 1pni
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_1pni_convection_box-reference.vtp
                               ${CMAKE_CURRENT_BINARY_DIR}/test_1pni_convection_box-00009.vtp
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pni_convection_box params.input -Problem.Name test_1pni_convection_box"
                       --zeroThreshold {"velocity_liq \(m/s\)":1e-15})

dumux_add_test(NAME test_1pni_convection_box_restart
              TARGET test_1pni_convection_box
              LABELS porousmediumflow 1p 1pni
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_1pni_convection_box-reference.vtp
                               ${CMAKE_CURRENT_BINARY_DIR}/test_1pni_convection_box_restart-00004.vtp
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pni_convection_box params.input -Problem.Name test_1pni_convection_box_restart -TimeLoop.DtInitial 1000 -Restart.Time 16768 -Restart.File test_1pni_convection_box-00006.vtp"
                       --zeroThreshold {"velocity_liq \(m/s\)":1e-15})

# the restart test has to run after the test that produces the corresponding vtu file
set_tests_properties(test_1pni_convection_box_restart PROPERTIES DEPENDS test_1pni_convection_box)

dumux_add_test(NAME test_1pni_convection_tpfa
              LABELS porousmediumflow 1p 1pni
              SOURCES ../main.cc
              COMPILE_DEFINITIONS TYPETAG=OnePNIConvectionCCTpfa
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_1pni_convection_cc-reference.vtp
                               ${CMAKE_CURRENT_BINARY_DIR}/test_1pni_convection_tpfa-00009.vtp
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pni_convection_tpfa params.input -Problem.Name test_1pni_convection_tpfa")

# the shared main.cc in the parent folder includes the properties.hh of this folder
foreach(target test_1pni_convection_box test_1pni_convection_tpfa)
  target_include_directories(${target} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
endforeach()
