# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links(FILES params.input)

dumux_add_test(NAME test_richardsnc_box
              LABELS porousmediumflow richardsnc
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=RichardsWellTracerBox
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_richardsnc_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_richardsnc_box-00035.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_richardsnc_box params.input -Problem.Name test_richardsnc_box")

dumux_add_test(NAME test_richardsnc_tpfa
              LABELS porousmediumflow richardsnc
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=RichardsWellTracerCC
              COMMAND ${DUMUX_RUNTEST}
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_richardsnc_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_richardsnc_tpfa-00033.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_richardsnc_tpfa params.input -Problem.Name test_richardsnc_tpfa")
