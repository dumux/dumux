# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()

dumux_add_test(NAME test_boussinesq_onesided_rb_box
              LABELS porousmediumflow boussinesq
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=BoussinesqOneSidedRBBox
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_boussinesq_onesided_rb_box
              CMD_ARGS params.input -Problem.Name test_boussinesq_onesided_rb_box)

dumux_add_test(NAME test_boussinesq_onesided_rb_box_convergence
              LABELS porousmediumflow boussinesq
              TARGET test_boussinesq_onesided_rb_box
              TIMEOUT 100000
              COMMAND ${Python3_EXECUTABLE}
              CMD_ARGS ${CMAKE_CURRENT_SOURCE_DIR}/convergence.py
                       ${CMAKE_CURRENT_BINARY_DIR}/test_boussinesq_onesided_rb_box
                       5)