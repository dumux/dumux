# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links(FILES params.input)

dumux_add_test(NAME test_3pwateroil_sagd_box
              LABELS porousmediumflow 3pwateroil
              SOURCES main.cc
              TIMEOUT 1500
              COMPILE_DEFINITIONS TYPETAG=ThreePWaterOilSagdBox
              COMMAND ${DUMUX_RUNTEST}
                CMD_ARGS --script fuzzy
            --files ${CMAKE_SOURCE_DIR}/test/references/test_3pwateroil_sagd_box-reference.vtu
                   ${CMAKE_CURRENT_BINARY_DIR}/test_3pwateroil_sagd_box-00001.vtu
            --command "${CMAKE_CURRENT_BINARY_DIR}/test_3pwateroil_sagd_box params.input -Problem.Name test_3pwateroil_sagd_box")
