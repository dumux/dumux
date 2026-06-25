# SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links()

dumux_add_test(NAME test_boussinesq_pressure_box
              LABELS porousmediumflow boussinesq
              SOURCES main.cc
              COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_boussinesq_pressure_box
              CMD_ARGS params.input -Problem.Name test_boussinesq_pressure_box)
